/*
 * SPDX-License-Identifier: GPL-2.0-only
 *
 * System-level evaluation (ns-3.48 + 5G-LENA v5.0) of handover interruption in
 * ultra-dense NR small-cell deployments, companion to the decision-delay model and the
 * outage-guarded A3 trigger (OG-A3).
 *
 *  - N x N square grid of gNBs wrapped on a torus (SquareWraparoundModel): no edge effects
 *  - n78 (3.5 GHz), numerology 1, 40 MHz, TDD DDDDU, real RRC, RLC UM, X2 handover
 *  - TR 38.901 UMi-StreetCanyon path loss with spatially correlated shadowing (ns-3) and a
 *    spatially consistent LOS state (OuLosConditionModel, exponential correlation, 50 m)
 *  - large-scale fading only (no small-scale fading), as in the trigger-level model
 *  - DL UDP CBR (1024 B / 5 ms) per UE, random-walk mobility (single or mixed speeds)
 *  - policies: A3-RSRP (Hys/TTT) and OG-A3 = A3 (Hys 3 dB, TTT 256 ms) + fast SINR guard
 *    + return block
 *
 * Per-HO records: Texec, Tpre = [tstart - tlast - Tpkt]+, Tpost, Tint, lost packets,
 * ping-pong (1 s / 3 s), trigger type and the measured decision delay TD = tstart - th,
 * where th is the instant the unfiltered RSRP difference first/last exceeded Hys.
 *
 * v3 additions (revision of the manuscript):
 *  - radio-level pre-HO gap Tpre,radio = [tstart - to]+, where to is the onset of the serving
 *    PDCCH-SINR outage (SINR < gammaMin) that is still ongoing at tstart, and the outage-margin
 *    time K = to - th of the HOs that suffer such a gap
 *  - MAC-level gap tstart - (last correctly decoded DL transport block from the source)
 *  - policy a3rb: A3 with the same return block as OG-A3 (HO back to the previous cell blocked
 *    for returnBlock seconds), implemented as a filter on the A3 measurement report
 *  - RLC UM options: t-Reordering value and out-of-order delivery
 */

#include "ns3/antenna-module.h"
#include "ns3/applications-module.h"
#include "ns3/core-module.h"
#include "ns3/flow-monitor-module.h"
#include "ns3/internet-module.h"
#include "ns3/mobility-module.h"
#include "ns3/network-module.h"
#include "ns3/nr-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/propagation-module.h"
#include "ns3/wraparound-model.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <vector>

using namespace ns3;

NS_LOG_COMPONENT_DEFINE("NrUdnHandoverV3");

// ---------------------------------------------------------------- square torus wraparound
class SquareWraparoundModel : public WraparoundModel
{
  public:
    static TypeId GetTypeId()
    {
        static TypeId tid = TypeId("ns3::SquareWraparoundModelV3")
                                .SetParent<WraparoundModel>()
                                .AddConstructor<SquareWraparoundModel>()
                                .AddAttribute("Side",
                                              "Side length of the square torus (m)",
                                              DoubleValue(300.0),
                                              MakeDoubleAccessor(&SquareWraparoundModel::m_side),
                                              MakeDoubleChecker<double>(0.0));
        return tid;
    }

    Vector3D GetVirtualPosition(const Vector3D tx, const Vector3D rx) const override
    {
        auto wrap = [this](double d) {
            double h = m_side / 2.0;
            double w = std::fmod(d + h, m_side);
            if (w < 0)
            {
                w += m_side;
            }
            return w - h;
        };
        return Vector3D(rx.x + wrap(tx.x - rx.x), rx.y + wrap(tx.y - rx.y), tx.z);
    }

  private:
    double m_side{300.0};
};

NS_OBJECT_ENSURE_REGISTERED(SquareWraparoundModel);

// ---------------------------------------------------------------- spatially consistent LOS
// LOS state of a link = 1{Phi(u) < PrLOS(d2D)}, u an exponentially correlated Gaussian process
// along the relative displacement of the link (decorrelation distance 50 m), as in the
// trigger-level Monte-Carlo simulator. TR 38.901 UMi-SC LOS probability.
class OuLosConditionModel : public ChannelConditionModel
{
  public:
    static TypeId GetTypeId()
    {
        static TypeId tid =
            TypeId("ns3::OuLosConditionModelV3")
                .SetParent<ChannelConditionModel>()
                .AddConstructor<OuLosConditionModel>()
                .AddAttribute("DecorrelationDistance",
                              "Decorrelation distance of the LOS process (m)",
                              DoubleValue(50.0),
                              MakeDoubleAccessor(&OuLosConditionModel::m_dcor),
                              MakeDoubleChecker<double>(0.0));
        return tid;
    }

    OuLosConditionModel()
    {
        m_norm = CreateObject<NormalRandomVariable>();
    }

    Ptr<ChannelCondition> GetChannelCondition(Ptr<const MobilityModel> a,
                                              Ptr<const MobilityModel> b) const override
    {
        uint32_t ia = a->GetObject<Node>()->GetId();
        uint32_t ib = b->GetObject<Node>()->GetId();
        uint64_t key = (static_cast<uint64_t>(std::min(ia, ib)) << 32) | std::max(ia, ib);
        Vector pa = a->GetPosition();
        Vector pb = b->GetPosition();
        // relative vector oriented from the lower to the higher node id
        Vector rel = ia < ib ? Vector(pb.x - pa.x, pb.y - pa.y, 0) : Vector(pa.x - pb.x, pa.y - pb.y, 0);
        auto it = m_state.find(key);
        double u;
        if (it == m_state.end())
        {
            u = m_norm->GetValue();
        }
        else
        {
            double disp = std::hypot(rel.x - it->second.rel.x, rel.y - it->second.rel.y);
            double r = std::exp(-disp / m_dcor);
            u = r * it->second.u + std::sqrt(1 - r * r) * m_norm->GetValue();
        }
        m_state[key] = {u, rel};
        double d2 = std::max(std::hypot(rel.x, rel.y), 1e-3);
        double pLos = d2 <= 18.0 ? 1.0 : 18.0 / d2 + std::exp(-d2 / 36.0) * (1 - 18.0 / d2);
        double phi = 0.5 * std::erfc(-u / std::sqrt(2.0));
        auto cond = phi < pLos ? ChannelCondition::LOS : ChannelCondition::NLOS;
        return CreateObject<ChannelCondition>(cond, ChannelCondition::O2O);
    }

    int64_t AssignStreams(int64_t stream) override
    {
        m_norm->SetStream(stream);
        return 1;
    }

  private:
    struct State
    {
        double u;
        Vector rel;
    };

    mutable std::map<uint64_t, State> m_state;
    Ptr<NormalRandomVariable> m_norm;
    double m_dcor{50.0};
};

NS_OBJECT_ENSURE_REGISTERED(OuLosConditionModel);

// ---------------------------------------------------------------- global state
namespace
{

struct HoRec
{
    uint32_t ue;
    double speedKmh;
    uint16_t src;
    uint16_t tgt;
    bool guard = false;
    double tStart;
    double tEnd = -1;
    double tLast = -1;  // last DL packet before HO start
    double tFirst = -1; // first DL packet after HO end
    int64_t seqLast = -1;
    int64_t lost = -1;
    bool pp1 = false;
    bool pp3 = false;
    double tdFirst = NAN; // tStart - first unfiltered A3 entry
    double tdLast = NAN;  // tStart - last unfiltered A3 entry
    double tReport = NAN; // last measurement report received by the source before tStart
    double sinrAtLast = NAN; // serving PDCCH SINR (dB) when the last packet before the HO arrived
    double sinrAtStart = NAN; // serving PDCCH SINR (dB) at HandoverStart
    double tOut = NAN;        // onset of the serving outage still ongoing at HandoverStart
    double tTbOk = NAN;       // last correctly decoded DL transport block from the source
    double kFirst = NAN;      // outage-margin time to - th (first unfiltered A3 entry)
    double kLast = NAN;       // outage-margin time to - th (last unfiltered A3 entry)
};

struct Ue
{
    Ptr<Node> node;
    Ptr<NrUeNetDevice> dev;
    double speedKmh = 0;
    // packets
    int64_t lastSeq = -1;
    double lastRx = -1;
    uint64_t rx = 0;
    int32_t pendingHo = -1; // index into g_hos awaiting first packet after HO
    bool inHo = false;
    // ping-pong
    uint16_t prevCell = 0;
    double prevHoTime = -1e9;
    // guard
    double sinrEmaDb = NAN;
    double sinrEmaT = -1;
    std::vector<double> rsrpEma; // per cell (dB), fast L1 estimate
    bool hoPending = false;
    double hoPendingSince = -1;
    bool guardIssued = false;
    // unfiltered A3 entry instants per cell (since last HO)
    std::vector<double> thFirst;
    std::vector<double> thLast;
    double lastReport = NAN;
    double lastCtrlSinrDb = NAN;
    double sinrAtLastRx = NAN;
    double outSince = NAN; // onset of the current serving outage (NaN if not in outage)
    double lastTbOk = NAN; // last correctly decoded DL transport block from the serving cell
};

std::vector<Ue> g_ue;
std::vector<HoRec> g_hos;
std::map<uint64_t, uint32_t> g_imsiToUe;
std::map<uint16_t, Ptr<NrGnbNetDevice>> g_cellToGnb;
std::map<uint16_t, uint32_t> g_cellToIdx;
std::vector<uint16_t> g_idxToCell;
NodeContainer g_gnbNodes;
Ptr<PropagationLossModel> g_pl;
Ptr<WraparoundModel> g_wrap;

// configuration used by callbacks
double g_warm = 2.0;
double g_hys = 3.0;
double g_pktInt = 0.005;
double g_gammaMin = -6.0;
bool g_guard = false;
double g_guardThr = -5.5;
double g_guardHys = 1.0;
double g_guardTavg = 0.02;
double g_guardSigma = 1.0;
double g_guardReportDelay = 0.005;
double g_returnBlock = 1.0;
bool g_trackTh = true;
double g_sampleDt = 0.005;
uint32_t g_hoFailures = 0;
uint32_t g_rlf = 0;
uint64_t g_outSamples = 0;
uint64_t g_totSamples = 0;
uint32_t g_guardFired = 0;
Ptr<NormalRandomVariable> g_noise;

void
UeHoStart(uint32_t u, uint64_t, uint16_t cellId, uint16_t, uint16_t target)
{
    double now = Simulator::Now().GetSeconds();
    Ue& s = g_ue[u];
    s.inHo = true;
    HoRec r;
    r.ue = u;
    r.speedKmh = s.speedKmh;
    r.src = cellId;
    r.tgt = target;
    r.tStart = now;
    r.tLast = s.lastRx;
    r.seqLast = s.lastSeq;
    r.guard = s.guardIssued;
    r.pp1 = (target == s.prevCell) && (now - s.prevHoTime <= 1.0);
    r.pp3 = (target == s.prevCell) && (now - s.prevHoTime <= 3.0);
    r.tReport = s.lastReport;
    r.sinrAtLast = s.sinrAtLastRx;
    r.sinrAtStart = s.lastCtrlSinrDb;
    r.tOut = s.outSince;
    r.tTbOk = s.lastTbOk;
    if (g_trackTh && g_cellToIdx.count(target))
    {
        uint32_t c = g_cellToIdx[target];
        r.tdFirst = now - s.thFirst[c];
        r.tdLast = now - s.thLast[c];
        r.kFirst = r.tOut - s.thFirst[c];
        r.kLast = r.tOut - s.thLast[c];
    }
    s.outSince = NAN;
    s.guardIssued = false;
    if (now >= g_warm)
    {
        g_hos.push_back(r);
        s.pendingHo = static_cast<int32_t>(g_hos.size() - 1);
    }
    else
    {
        s.pendingHo = -1;
    }
}

void
UeHoEndOk(uint32_t u, uint64_t, uint16_t cellId, uint16_t)
{
    double now = Simulator::Now().GetSeconds();
    Ue& s = g_ue[u];
    s.inHo = false;
    s.hoPending = false;
    if (s.pendingHo >= 0)
    {
        HoRec& r = g_hos[s.pendingHo];
        r.tEnd = now;
        s.prevCell = r.src;
    }
    else
    {
        s.prevCell = 0;
    }
    s.prevHoTime = now;
    std::fill(s.thFirst.begin(), s.thFirst.end(), std::numeric_limits<double>::infinity());
    std::fill(s.thLast.begin(), s.thLast.end(), std::numeric_limits<double>::infinity());
}

void
UeHoEndError(uint32_t u, uint64_t, uint16_t, uint16_t)
{
    Ue& s = g_ue[u];
    s.inHo = false;
    s.hoPending = false;
    if (Simulator::Now().GetSeconds() >= g_warm)
    {
        ++g_hoFailures;
    }
    s.pendingHo = -1;
}

void
UeRlf(uint32_t u, uint64_t, uint16_t, uint16_t)
{
    g_ue[u].hoPending = false;
    if (Simulator::Now().GetSeconds() >= g_warm)
    {
        ++g_rlf;
    }
}

void
GnbMeasReport(uint64_t imsi, uint16_t, uint16_t, NrRrcSap::MeasurementReport)
{
    auto it = g_imsiToUe.find(imsi);
    if (it != g_imsiToUe.end())
    {
        g_ue[it->second].lastReport = Simulator::Now().GetSeconds();
    }
}

void
DlCtrlSinr(uint32_t u, uint16_t cellId, uint16_t, double sinrLin, uint16_t)
{
    Ue& s = g_ue[u];
    if (cellId != s.dev->GetRrc()->GetCellId())
    {
        return;
    }
    double now = Simulator::Now().GetSeconds();
    double db = 10 * std::log10(std::max(sinrLin, 1e-12));
    s.lastCtrlSinrDb = db;
    if (!g_guard)
    {
        return;
    }
    double est = db + g_guardSigma * g_noise->GetValue();
    if (std::isnan(s.sinrEmaDb))
    {
        s.sinrEmaDb = est;
    }
    else
    {
        double w = 1 - std::exp(-(now - s.sinrEmaT) / g_guardTavg);
        s.sinrEmaDb = (1 - w) * s.sinrEmaDb + w * est;
    }
    s.sinrEmaT = now;
}

void
UeRxTb(uint32_t u, RxPacketTraceParams p)
{
    Ue& s = g_ue[u];
    if (!p.m_corrupt && p.m_cellId == s.dev->GetRrc()->GetCellId())
    {
        s.lastTbOk = Simulator::Now().GetSeconds();
    }
}

void
DlPacketRx(uint32_t u, Ptr<const Packet> p, const Address&)
{
    SeqTsHeader h;
    Ptr<Packet> c = p->Copy();
    c->RemoveHeader(h);
    int64_t seq = h.GetSeq();
    double now = Simulator::Now().GetSeconds();
    Ue& s = g_ue[u];
    if (s.pendingHo >= 0 && g_hos[s.pendingHo].tEnd >= 0 && !s.inHo)
    {
        HoRec& r = g_hos[s.pendingHo];
        r.tFirst = now;
        r.lost = std::max<int64_t>(0, seq - r.seqLast - 1);
        s.pendingHo = -1;
    }
    s.lastSeq = std::max(s.lastSeq, seq);
    s.lastRx = now;
    s.sinrAtLastRx = s.lastCtrlSinrDb;
    ++s.rx;
}

void
IssueGuardHo(uint32_t u, uint16_t servCell, uint16_t target)
{
    Ue& s = g_ue[u];
    Ptr<NrUeRrc> rrc = s.dev->GetRrc();
    if (rrc->GetCellId() != servCell || rrc->GetState() != NrUeRrc::CONNECTED_NORMALLY)
    {
        s.hoPending = false;
        return;
    }
    Ptr<NrGnbRrc> grrc = g_cellToGnb[servCell]->GetRrc();
    uint16_t rnti = rrc->GetRnti();
    if (!grrc->HasUeManager(rnti) ||
        grrc->GetUeManager(rnti)->GetState() != NrUeManager::CONNECTED_NORMALLY)
    {
        s.hoPending = false;
        return;
    }
    s.guardIssued = true;
    ++g_guardFired;
    grrc->SendHandoverRequest(rnti, target);
}

void
Turn(Ptr<ConstantVelocityMobilityModel> m,
     double v,
     double isd,
     Ptr<UniformRandomVariable> uni,
     double firstLeg)
{
    double a = uni->GetValue(0, 2 * M_PI);
    m->SetVelocity(Vector(v * std::cos(a), v * std::sin(a), 0));
    if (v > 0)
    {
        double leg = firstLeg > 0 ? firstLeg : isd;
        Simulator::Schedule(Seconds(leg / v), &Turn, m, v, isd, uni, 0.0);
    }
}

// 5 ms loop: outage sampling, unfiltered A3-entry tracking, fast L1 RSRP and the guard
void
Sample()
{
    double now = Simulator::Now().GetSeconds();
    uint32_t nCells = g_gnbNodes.GetN();
    for (uint32_t u = 0; u < g_ue.size(); ++u)
    {
        Ue& s = g_ue[u];
        Ptr<NrUeRrc> rrc = s.dev->GetRrc();
        uint16_t serv = rrc->GetCellId();
        if (now >= g_warm)
        {
            ++g_totSamples;
            bool bad = s.inHo || std::isnan(s.lastCtrlSinrDb) || s.lastCtrlSinrDb < g_gammaMin;
            g_outSamples += bad ? 1 : 0;
        }
        // onset of the current radio outage of the serving link (radio-level pre-HO gap)
        bool radioBad = !s.inHo && !std::isnan(s.lastCtrlSinrDb) && s.lastCtrlSinrDb < g_gammaMin;
        if (!radioBad)
        {
            s.outSince = NAN;
        }
        else if (std::isnan(s.outSince))
        {
            s.outSince = now;
        }
        if (!g_cellToIdx.count(serv) || (!g_trackTh && !g_guard))
        {
            continue;
        }
        uint32_t sc = g_cellToIdx[serv];
        Ptr<MobilityModel> um = s.node->GetObject<MobilityModel>();
        std::vector<double> rsrp(nCells);
        for (uint32_t c = 0; c < nCells; ++c)
        {
            Ptr<MobilityModel> gm = g_gnbNodes.Get(c)->GetObject<MobilityModel>();
            Ptr<MobilityModel> vm = g_wrap->GetVirtualMobilityModel(gm, um);
            rsrp[c] = g_pl->CalcRxPower(0.0, vm, um);
        }
        if (g_trackTh && !s.inHo)
        {
            for (uint32_t c = 0; c < nCells; ++c)
            {
                bool cond = (c != sc) && (rsrp[c] - rsrp[sc] > g_hys);
                if (cond)
                {
                    s.thFirst[c] = std::min(s.thFirst[c], now);
                    if (std::isinf(s.thLast[c]))
                    {
                        s.thLast[c] = now;
                    }
                }
                else
                {
                    s.thLast[c] = std::numeric_limits<double>::infinity();
                }
            }
        }
        if (!g_guard)
        {
            continue;
        }
        double w = 1 - std::exp(-g_sampleDt / g_guardTavg);
        for (uint32_t c = 0; c < nCells; ++c)
        {
            double meas = rsrp[c] + g_guardSigma * g_noise->GetValue();
            s.rsrpEma[c] = std::isnan(s.rsrpEma[c]) ? meas : (1 - w) * s.rsrpEma[c] + w * meas;
        }
        if (s.hoPending && now - s.hoPendingSince > 1.0)
        {
            s.hoPending = false; // request was ignored or failed silently
        }
        if (s.inHo || s.hoPending || std::isnan(s.sinrEmaDb) || s.sinrEmaDb >= g_guardThr ||
            rrc->GetState() != NrUeRrc::CONNECTED_NORMALLY)
        {
            continue;
        }
        int32_t best = -1;
        double bestV = -1e9;
        for (uint32_t c = 0; c < nCells; ++c)
        {
            uint16_t cid = g_idxToCell[c];
            if (c == sc || (cid == s.prevCell && now - s.prevHoTime <= g_returnBlock))
            {
                continue;
            }
            if (s.rsrpEma[c] > bestV)
            {
                bestV = s.rsrpEma[c];
                best = c;
            }
        }
        if (best >= 0 && bestV > s.rsrpEma[sc] + g_guardHys)
        {
            uint16_t tcid = g_idxToCell[best];
            s.hoPending = true;
            s.hoPendingSince = now;
            Simulator::Schedule(Seconds(g_guardReportDelay), &IssueGuardHo, u, serv, tcid);
        }
    }
    Simulator::Schedule(Seconds(g_sampleDt), &Sample);
}

} // namespace

// ---------------------------------------------------------------- A3 + return block
// Standard A3-RSRP algorithm whose measurement reports are filtered so that the UE cannot be
// handed back to the cell it left less than g_returnBlock seconds ago (same rule as OG-A3).
class NrA3ReturnBlockHandoverAlgorithm : public NrA3RsrpHandoverAlgorithm
{
  public:
    static TypeId GetTypeId()
    {
        static TypeId tid = TypeId("ns3::NrA3ReturnBlockHandoverAlgorithmV3")
                                .SetParent<NrA3RsrpHandoverAlgorithm>()
                                .AddConstructor<NrA3ReturnBlockHandoverAlgorithm>();
        return tid;
    }

    static std::map<const NrHandoverAlgorithm*, uint16_t>& CellOf()
    {
        static std::map<const NrHandoverAlgorithm*, uint16_t> m;
        return m;
    }

  protected:
    void DoReportUeMeas(uint16_t rnti, NrRrcSap::MeasResults mr) override
    {
        auto it = CellOf().find(this);
        if (it != CellOf().end())
        {
            double now = Simulator::Now().GetSeconds();
            for (const auto& s : g_ue)
            {
                Ptr<NrUeRrc> rrc = s.dev->GetRrc();
                if (rrc->GetCellId() == it->second && rrc->GetRnti() == rnti)
                {
                    if (s.prevCell != 0 && now - s.prevHoTime <= g_returnBlock)
                    {
                        mr.measResultListEutra.remove_if(
                            [&s](const NrRrcSap::MeasResultEutra& e) {
                                return e.physCellId == s.prevCell;
                            });
                        mr.haveMeasResultNeighCells = !mr.measResultListEutra.empty();
                    }
                    break;
                }
            }
        }
        NrA3RsrpHandoverAlgorithm::DoReportUeMeas(rnti, mr);
    }
};

NS_OBJECT_ENSURE_REGISTERED(NrA3ReturnBlockHandoverAlgorithm);

int
main(int argc, char* argv[])
{
    uint16_t gridN = 6;
    double isd = 50.0;
    uint32_t uesPerCell = 2;
    std::string speeds = "30"; // comma-separated km/h; UEs split equally (mixed population)
    double simTime = 30.0;
    std::string policy = "a3"; // a3 | og
    double ttt = 100;          // ms (A3); OG-A3 uses ogTtt
    double ogTtt = 256;
    double tm = 40; // ms, UE measurement (L1) period
    double txPower = 30;
    double bw = 40e6;
    double fc = 3.5e9;
    uint16_t numerology = 1;
    std::string pattern = "DL|DL|DL|DL|UL|";
    std::string prefix = "nrudn";
    uint32_t run = 1;
    bool rlcOoo = false;
    double rlcReorderMs = 100;

    CommandLine cmd(__FILE__);
    cmd.AddValue("gridN", "gNB grid size N (N x N on a torus)", gridN);
    cmd.AddValue("isd", "Inter-site distance (m)", isd);
    cmd.AddValue("uesPerCell", "UEs per cell", uesPerCell);
    cmd.AddValue("speeds", "Comma-separated UE speeds (km/h), split equally", speeds);
    cmd.AddValue("simTime", "Simulation time (s)", simTime);
    cmd.AddValue("warm", "Warm-up time excluded from KPIs (s)", g_warm);
    cmd.AddValue("policy", "a3 | a3rb (A3 + return block) | og", policy);
    cmd.AddValue("rlcOoo", "RLC UM out-of-order delivery", rlcOoo);
    cmd.AddValue("rlcReorderMs", "RLC UM t-Reordering (ms)", rlcReorderMs);
    cmd.AddValue("hys", "A3 hysteresis (dB)", g_hys);
    cmd.AddValue("ttt", "A3 time-to-trigger (ms)", ttt);
    cmd.AddValue("ogTtt", "OG-A3 robust A3 time-to-trigger (ms)", ogTtt);
    cmd.AddValue("tm", "UE L1 measurement period (ms)", tm);
    cmd.AddValue("guardThr", "OG-A3 guard SINR threshold (dB)", g_guardThr);
    cmd.AddValue("guardHys", "OG-A3 guard neighbour offset (dB)", g_guardHys);
    cmd.AddValue("guardTavg", "OG-A3 guard averaging window (s)", g_guardTavg);
    cmd.AddValue("guardSigma", "OG-A3 guard estimation noise (dB)", g_guardSigma);
    cmd.AddValue("guardDelay", "OG-A3 L1 report delay (s)", g_guardReportDelay);
    cmd.AddValue("returnBlock", "OG-A3 return block (s)", g_returnBlock);
    cmd.AddValue("trackTh", "Track unfiltered A3 entry for model validation", g_trackTh);
    cmd.AddValue("gammaMin", "Outage SINR threshold (dB)", g_gammaMin);
    cmd.AddValue("txPower", "gNB Tx power (dBm)", txPower);
    cmd.AddValue("pktInterval", "DL UDP inter-packet time (s)", g_pktInt);
    cmd.AddValue("prefix", "Output prefix", prefix);
    cmd.AddValue("run", "RNG run", run);
    cmd.Parse(argc, argv);

    RngSeedManager::SetSeed(1);
    RngSeedManager::SetRun(run);
    g_guard = (policy == "og");
    g_noise = CreateObject<NormalRandomVariable>();

    std::vector<double> speedList;
    {
        std::stringstream ss(speeds);
        std::string tok;
        while (std::getline(ss, tok, ','))
        {
            speedList.push_back(std::stod(tok));
        }
    }

    const uint32_t nCells = gridN * gridN;
    const uint32_t nUes = nCells * uesPerCell;
    const double side = gridN * isd;

    Config::SetDefault("ns3::NrUePhy::UeMeasurementsFilterPeriod", TimeValue(MilliSeconds(tm)));
    Config::SetDefault("ns3::NrGnbRrc::QosFlowToRlcMapping", StringValue("RlcUmAlways"));
    Config::SetDefault("ns3::NrRlcUm::MaxTxBufferSize", UintegerValue(1024 * 1024));
    Config::SetDefault("ns3::NrRlcUm::OutOfOrderDelivery", BooleanValue(rlcOoo));
    Config::SetDefault("ns3::NrRlcUm::ReorderingTimer", TimeValue(MilliSeconds(rlcReorderMs)));
    Config::SetDefault("ns3::UdpClient::Interval", TimeValue(Seconds(g_pktInt)));
    Config::SetDefault("ns3::UdpClient::MaxPackets", UintegerValue(0));
    Config::SetDefault("ns3::UdpClient::PacketSize", UintegerValue(1024));

    Ptr<NrHelper> nrHelper = CreateObject<NrHelper>();
    Ptr<NrPointToPointEpcHelper> epcHelper = CreateObject<NrPointToPointEpcHelper>();
    nrHelper->SetEpcHelper(epcHelper);
    nrHelper->SetAttribute("UseIdealRrc", BooleanValue(false));
    epcHelper->SetAttribute("S1uLinkDelay", TimeValue(MilliSeconds(2)));
    nrHelper->SetHandoverAlgorithmType(policy == "a3rb" ? "ns3::NrA3ReturnBlockHandoverAlgorithmV3"
                                                        : "ns3::NrA3RsrpHandoverAlgorithm");
    nrHelper->SetHandoverAlgorithmAttribute("Hysteresis", DoubleValue(g_hys));
    nrHelper->SetHandoverAlgorithmAttribute(
        "TimeToTrigger",
        TimeValue(MilliSeconds(policy == "og" ? ogTtt : ttt)));
    // SRS disabled: DL-only traffic, CQI from PDCCH; avoids a 5G-LENA v5.0 SRS-collision abort
    // that occurs after frequent handovers ("Not allowed state for starting SRS reception")
    nrHelper->SetSchedulerAttribute("EnableSrsInUlSlots", BooleanValue(false));
    nrHelper->SetSchedulerAttribute("EnableSrsInFSlots", BooleanValue(false));
    nrHelper->SetUeAntennaTypeId(IsotropicAntennaModel::GetTypeId().GetName());
    nrHelper->SetGnbAntennaTypeId(IsotropicAntennaModel::GetTypeId().GetName());
    nrHelper->SetGnbPhyAttribute("TxPower", DoubleValue(txPower));
    nrHelper->SetGnbPhyAttribute("Numerology", UintegerValue(numerology));
    nrHelper->SetGnbPhyAttribute("Pattern", StringValue(pattern));
    nrHelper->SetUePhyAttribute("TxPower", DoubleValue(23));
    nrHelper->SetUePhyAttribute("NoiseFigure", DoubleValue(9));

    // channel: UMi-SC large-scale only, torus wraparound, spatially consistent LOS
    auto wrap = CreateObject<SquareWraparoundModel>();
    wrap->SetAttribute("Side", DoubleValue(side));
    g_wrap = wrap;
    Ptr<NrChannelHelper> channelHelper = CreateObject<NrChannelHelper>();
    channelHelper->ConfigureFactories("UMi", "Default", "ThreeGpp");
    channelHelper->SetPathlossAttribute("ShadowingEnabled", BooleanValue(true));
    channelHelper->SetWraparoundModel(wrap);
    CcBwpCreator ccBwpCreator;
    CcBwpCreator::SimpleOperationBandConf bandConf(fc, bw, 1);
    OperationBandInfo band = ccBwpCreator.CreateOperationBandContiguousCc(bandConf);
    channelHelper->AssignChannelsToBands({band}, NrChannelHelper::INIT_PROPAGATION);
    BandwidthPartInfoPtrVector allBwps = CcBwpCreator::GetAllBwps({band});
    auto losModel = CreateObject<OuLosConditionModel>();
    g_pl = band.m_cc[0]->m_bwp[0]->GetChannel()->GetPropagationLossModel();
    g_pl->SetAttribute("ChannelConditionModel", PointerValue(losModel));

    // nodes
    g_gnbNodes.Create(nCells);
    NodeContainer ueNodes;
    ueNodes.Create(nUes);
    Ptr<ListPositionAllocator> gpos = CreateObject<ListPositionAllocator>();
    for (uint32_t i = 0; i < nCells; ++i)
    {
        gpos->Add(Vector((i % gridN + 0.5) * isd, (i / gridN + 0.5) * isd, 10.0));
    }
    MobilityHelper gmob;
    gmob.SetMobilityModel("ns3::ConstantPositionMobilityModel");
    gmob.SetPositionAllocator(gpos);
    gmob.Install(g_gnbNodes);

    // UE mobility: random walk on the torus (constant speed, new uniformly random direction
    // every isd metres). Positions are not bounded; the wraparound model maps every link to
    // its closest virtual image, so there are no edge effects and no boundary reflections.
    g_ue.resize(nUes);
    Ptr<UniformRandomVariable> uni = CreateObject<UniformRandomVariable>();
    MobilityHelper um;
    um.SetMobilityModel("ns3::ConstantVelocityMobilityModel");
    um.Install(ueNodes);
    for (uint32_t u = 0; u < nUes; ++u)
    {
        double v = speedList[u % speedList.size()];
        g_ue[u].speedKmh = v;
        ueNodes.Get(u)->GetObject<MobilityModel>()->SetPosition(
            Vector(uni->GetValue(0, side), uni->GetValue(0, side), 1.5));
        Turn(ueNodes.Get(u)->GetObject<ConstantVelocityMobilityModel>(), v / 3.6, isd, uni,
             uni->GetValue(0, isd));
    }

    NetDeviceContainer gnbDevs = nrHelper->InstallGnbDevice(g_gnbNodes, allBwps);
    NetDeviceContainer ueDevs = nrHelper->InstallUeDevice(ueNodes, allBwps);
    int64_t stream = 1;
    stream += nrHelper->AssignStreams(gnbDevs, stream);
    stream += nrHelper->AssignStreams(ueDevs, stream);
    stream += losModel->AssignStreams(stream);
    g_noise->SetStream(stream++);

    auto [remoteHost, remoteHostAddr] = epcHelper->SetupRemoteHost("100Gb/s", 1500, Seconds(0));
    InternetStackHelper internet;
    internet.Install(ueNodes);
    Ipv4InterfaceContainer ueIps = epcHelper->AssignUeIpv4Address(ueDevs);
    nrHelper->AttachToClosestGnb(ueDevs, gnbDevs);
    nrHelper->AddX2Interface(g_gnbNodes);

    for (uint32_t i = 0; i < gnbDevs.GetN(); ++i)
    {
        Ptr<NrGnbNetDevice> g = gnbDevs.Get(i)->GetObject<NrGnbNetDevice>();
        g_cellToGnb[g->GetCellId()] = g;
        g_cellToIdx[g->GetCellId()] = i;
        g_idxToCell.push_back(g->GetCellId());
        g->GetRrc()->TraceConnectWithoutContext("RecvMeasurementReport",
                                                MakeCallback(&GnbMeasReport));
        PointerValue hoAlgo;
        g->GetAttribute("NrHandoverAlgorithm", hoAlgo);
        NrA3ReturnBlockHandoverAlgorithm::CellOf()[PeekPointer(hoAlgo.Get<NrHandoverAlgorithm>())] =
            g->GetCellId();
    }

    // applications
    const uint16_t port = 10000;
    ApplicationContainer clients;
    ApplicationContainer sinks;
    for (uint32_t u = 0; u < nUes; ++u)
    {
        Ue& s = g_ue[u];
        s.node = ueNodes.Get(u);
        s.dev = ueDevs.Get(u)->GetObject<NrUeNetDevice>();
        s.thFirst.assign(nCells, std::numeric_limits<double>::infinity());
        s.thLast.assign(nCells, std::numeric_limits<double>::infinity());
        s.rsrpEma.assign(nCells, NAN);
        g_imsiToUe[s.dev->GetImsi()] = u;
        Ptr<NrUeRrc> rrc = s.dev->GetRrc();
        rrc->TraceConnectWithoutContext("HandoverStart", MakeBoundCallback(&UeHoStart, u));
        rrc->TraceConnectWithoutContext("HandoverEndOk", MakeBoundCallback(&UeHoEndOk, u));
        rrc->TraceConnectWithoutContext("HandoverEndError", MakeBoundCallback(&UeHoEndError, u));
        rrc->TraceConnectWithoutContext("RadioLinkFailure", MakeBoundCallback(&UeRlf, u));
        s.dev->GetPhy(0)->TraceConnectWithoutContext("DlCtrlSinr", MakeBoundCallback(&DlCtrlSinr, u));
        s.dev->GetPhy(0)->GetSpectrumPhy()->TraceConnectWithoutContext(
            "RxPacketTraceUe",
            MakeBoundCallback(&UeRxTb, u));

        Ptr<Ipv4StaticRouting> rt = Ipv4StaticRoutingHelper().GetStaticRouting(s.node->GetObject<Ipv4>());
        rt->SetDefaultRoute(epcHelper->GetUeDefaultGatewayAddress(), 1);
        UdpClientHelper c(ueIps.GetAddress(u), port);
        clients.Add(c.Install(remoteHost));
        PacketSinkHelper sk("ns3::UdpSocketFactory", InetSocketAddress(Ipv4Address::GetAny(), port));
        ApplicationContainer a = sk.Install(s.node);
        a.Get(0)->TraceConnectWithoutContext("Rx", MakeBoundCallback(&DlPacketRx, u));
        sinks.Add(a);
    }
    sinks.Start(Seconds(0.4));
    clients.Start(Seconds(0.5));
    clients.Stop(Seconds(simTime - 0.2));

    FlowMonitorHelper fmh;
    Ptr<FlowMonitor> fm = fmh.Install(ueNodes);
    fmh.Install(remoteHost);

    Simulator::Schedule(Seconds(0.5), &Sample);
    Simulator::Stop(Seconds(simTime));
    Simulator::Run();

    // ---------------------------------------------------------------- results
    fm->CheckForLostPackets();
    auto cls = DynamicCast<Ipv4FlowClassifier>(fmh.GetClassifier());
    uint64_t txP = 0, rxP = 0, rxB = 0;
    double dSum = 0;
    for (const auto& [id, st] : fm->GetFlowStats())
    {
        if (cls->FindFlow(id).destinationPort != port)
        {
            continue;
        }
        txP += st.txPackets;
        rxP += st.rxPackets;
        rxB += st.rxBytes;
        dSum += st.delaySum.GetSeconds();
    }
    const double active = simTime - 0.2 - 0.5;

    std::ofstream hf(prefix + "-handovers.csv");
    hf << "ue,speed_kmh,src,tgt,trigger,tStart,tExec_ms,tPre_ms,tPost_ms,tInt_ms,lost,pp1,pp3,"
          "td_first_ms,td_last_ms,tPrep_ms,sinrLast_db,sinrStart_db,tPreRadio_ms,k_first_ms,"
          "k_last_ms,macGap_ms\n";
    std::vector<double> exec, pre, post, inte, lost, preRadio, macGap;
    uint32_t pp1 = 0, pp3 = 0, nGuard = 0;
    for (const auto& r : g_hos)
    {
        double te = r.tEnd >= 0 ? (r.tEnd - r.tStart) * 1e3 : NAN;
        double tp = r.tLast >= 0 ? std::max(0.0, (r.tStart - r.tLast - g_pktInt) * 1e3) : NAN;
        double tpo = (r.tFirst >= 0 && r.tEnd >= 0) ? (r.tFirst - r.tEnd) * 1e3 : NAN;
        double ti = (r.tFirst >= 0 && r.tLast >= 0) ? (r.tFirst - r.tLast) * 1e3 : NAN;
        double tprep = std::isnan(r.tReport) ? NAN : (r.tStart - r.tReport) * 1e3;
        hf << r.ue << "," << r.speedKmh << "," << r.src << "," << r.tgt << ","
           << (r.guard ? "guard" : "a3") << "," << r.tStart << "," << te << "," << tp << ","
           << tpo << "," << ti << "," << r.lost << "," << r.pp1 << "," << r.pp3 << ","
           << r.tdFirst * 1e3 << "," << r.tdLast * 1e3 << "," << tprep << "," << r.sinrAtLast << ","
           << r.sinrAtStart << ",";
        double tpr = std::isnan(r.tOut) ? 0.0 : (r.tStart - r.tOut) * 1e3;
        double mg = std::isnan(r.tTbOk) ? NAN : (r.tStart - r.tTbOk) * 1e3;
        hf << tpr << "," << r.kFirst * 1e3 << "," << r.kLast * 1e3 << "," << mg << "\n";
        preRadio.push_back(tpr);
        if (!std::isnan(mg))
        {
            macGap.push_back(mg);
        }
        pp1 += r.pp1;
        pp3 += r.pp3;
        nGuard += r.guard;
        if (!std::isnan(te))
        {
            exec.push_back(te);
        }
        if (!std::isnan(tp))
        {
            pre.push_back(tp);
        }
        if (!std::isnan(tpo))
        {
            post.push_back(tpo);
        }
        if (!std::isnan(ti))
        {
            inte.push_back(ti);
            lost.push_back(r.lost);
        }
    }
    auto mean = [](const std::vector<double>& v) {
        return v.empty() ? NAN : std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    };
    double lostHo = std::accumulate(lost.begin(), lost.end(), 0.0);
    double n = g_hos.size();
    std::ofstream sf(prefix + "-summary.csv");
    sf << "run,isd,gridN,ues,speeds,policy,hys,ttt,tm,numHO,hoRate,texec_ms,tpre_ms,tpost_ms,"
          "tint_ms,lostPerHo,plr_pct,hoPlr_pct,pp1_pct,pp3_pct,guard_pct,hoFail,rlf,outage_pct,"
          "thr_Mbps,delay_ms,tpreRadio_ms,pRadioGap,macGap_ms\n";
    sf << run << "," << isd << "," << gridN << "," << nUes << ",\"" << speeds << "\"," << policy
       << "," << g_hys << "," << (policy == "og" ? ogTtt : ttt) << "," << tm << "," << n << ","
       << n / (nUes * (simTime - g_warm)) << "," << mean(exec) << "," << mean(pre) << ","
       << mean(post) << "," << mean(inte) << "," << mean(lost) << ","
       << (txP ? 100.0 * (txP - rxP) / txP : NAN) << "," << (txP ? 100.0 * lostHo / txP : NAN)
       << "," << (n ? 100.0 * pp1 / n : NAN) << "," << (n ? 100.0 * pp3 / n : NAN) << ","
       << (n ? 100.0 * nGuard / n : 0) << "," << g_hoFailures << "," << g_rlf << ","
       << (g_totSamples ? 100.0 * g_outSamples / g_totSamples : NAN) << ","
       << rxB * 8.0 / active / 1e6 / nUes << "," << (rxP ? dSum / rxP * 1e3 : NAN) << ","
       << mean(preRadio) << ","
       << (preRadio.empty() ? NAN
                            : std::count_if(preRadio.begin(),
                                            preRadio.end(),
                                            [](double x) { return x > 0; }) /
                                  double(preRadio.size()))
       << "," << mean(macGap) << "\n";

    std::cout << "HOs=" << n << " rate=" << n / (nUes * (simTime - g_warm)) << " Texec="
              << mean(exec) << " Tpre=" << mean(pre) << " Tint=" << mean(inte)
              << " PLR=" << (txP ? 100.0 * (txP - rxP) / txP : NAN) << "% PP1="
              << (n ? 100.0 * pp1 / n : NAN) << "% guard=" << (n ? 100.0 * nGuard / n : 0)
              << "% HOF=" << g_hoFailures << " RLF=" << g_rlf
              << " out=" << (g_totSamples ? 100.0 * g_outSamples / g_totSamples : NAN) << "%\n";
    Simulator::Destroy();
    return 0;
}
