#ifndef SND_PLANETOOLS_H
#define SND_PLANETOOLS_H        

#include <vector>
#include <algorithm>
#include <numeric>
#include <type_traits>

#include "TClonesArray.h"
#include "Scifi.h"
#include "MuFilter.h"
#include "sndConfiguration.h"
#include "sndVetoPlane.h"
#include "sndScifiPlane.h"
#include "sndUSPlane.h"
#include "sndDSPlane.h"

namespace snd {
    namespace analysis_tools {
        // Produce veto, scifi, us and ds planes from data 
        std::vector<VetoPlane> FillVeto(const Configuration &configuration, TClonesArray *mufi_hits, MuFilter *mufilter_geometry);
        std::vector<ScifiPlane> FillScifi(const Configuration &configuration, TClonesArray *sf_hits, Scifi *scifi_geometry);
        std::vector<USPlane> FillUS(const Configuration &configuration, TClonesArray *mufi_hits, MuFilter *mufilter_geometry, bool use_small_sipms=false);
        std::vector<DSPlane> FillDS(const Configuration &configuration, TClonesArray *mufi_hits, MuFilter *mufilter_geometry);
    
        struct Cluster {
            ROOT::Math::XYZPoint center;
            double radius;   // in 2D: circle radius, in 1D: half-length
            bool is2D;
        };

        template <typename H>
        std::vector<Cluster> ClustersPositions(std::vector<H> hits, double min_radius_x, double min_radius_y, double max_gap = 1.0) {

        static_assert(std::is_convertible_v<decltype(std::declval<const H>().HitPosition()), ROOT::Math::XYZPoint>,
                        "Hit class must have HitPosition()");

            std::vector<snd::analysis_tools::Cluster> clusters;
            if (hits.empty()) return clusters;

            // Check if any hit has both x and y valid (2D mixed hits)
            bool any2D = std::any_of(hits.begin(), hits.end(), [](const auto& h) {
                return !std::isnan(h.x) && !std::isnan(h.y);
            });

            // Lambda that sorts a group and clusters it along the given axis
            auto clusterAlong = [&](std::vector<H>& group, bool useX) {

                if (group.empty()) return;

                std::sort(group.begin(), group.end(), [useX](const auto& a, const auto& b) {
                    float aPos = useX ? a.x : a.y;
                    float bPos = useX ? b.x : b.y;
                    if (std::isnan(aPos)) return false;
                    if (std::isnan(bPos)) return true;
                    return aPos < bPos;
                });

                // std::cout << "DEBUG: sorted order:" << std::endl;
                // for (const auto& h :group) h.Print();

                size_t start{0};
                size_t N = group.size();
                while (start < N) {
                    size_t end = start;

                    while (end + 1 < N) {
                        float current_pos = useX ? group[end].x     : group[end].y;
                        float next_pos    = useX ? group[end + 1].x : group[end + 1].y;

                        if (std::abs(next_pos - current_pos) > max_gap) {
                                // std::cout << "DEBUG gap check: current_pos=" << current_pos 
                                // << " next_pos=" << next_pos 
                                // << " diff=" << std::abs(next_pos - current_pos)
                                // << " max_gap=" << max_gap << std::endl;
                                break;
                                }
    
                        ++end;
                    }

                    double sumx = 0.0, sumy = 0.0, sumz = 0.0;
                    int countX = 0, countY = 0, countZ = 0;
                    double min_pos = 1e9, max_pos = -1e9;

                    for (size_t i = start; i <= end; ++i) {
                        auto p = group[i].HitPosition();

                        if (!std::isnan(p.X())) {
                            sumx += p.X();
                            countX++;
                        }
                        if (!std::isnan(p.Y())) {
                            sumy += p.Y();
                            countY++;
                        }
                        if (!std::isnan(p.Z())) {
                            sumz += p.Z();
                            countZ++;
                        }

                        float pos = useX ? p.X() : p.Y();
                        if (!std::isnan(pos)) {
                            min_pos = std::min(min_pos, (double)pos);
                            max_pos = std::max(max_pos, (double)pos);
                        }
                    }

                    bool hasX = (countX > 0);
                    bool hasY = (countY > 0);
                    bool is2D = hasX && hasY;

                    ROOT::Math::XYZPoint center(
                        hasX ? sumx / countX : NAN,
                        hasY ? sumy / countY : NAN,
                        countZ > 0 ? sumz / countZ : NAN
                    );

                    double min_radius = useX ? min_radius_x : min_radius_y;
                    double radius = 0.0;

                    if (is2D) {
                        double max_dist_sq = 0.0;
                        for (size_t i = start; i <= end; ++i) {
                            auto p = group[i].HitPosition();
                            double dx = p.X() - center.X();
                            double dy = p.Y() - center.Y();
                            max_dist_sq = std::max(max_dist_sq, dx*dx + dy*dy);
                        }
                        radius = std::sqrt(max_dist_sq);
                        if (radius < min_radius) radius = min_radius;
                    } else if (hasX || hasY) {
                        radius = (max_pos - min_pos) / 2.0;
                        if (radius < min_radius) radius = min_radius;
                    }

                    // std::cout << "DEBUG: New Cluster from hits " << start << " to " << end
                    //         << " (Total hits in cluster: " << (end - start + 1) << ")" << std::endl;
                    // for (size_t i = start; i <= end; ++i) group[i].Print();

                    clusters.push_back({center, radius, is2D});
                    start = end + 1;
                }
            };

            if (any2D) {
                // All hits have y, some also have x — cluster together along y
                clusterAlong(hits, false);
            } else {
                // Hits are strictly x-only or y-only — split and cluster independently
                std::vector<H> xHits, yHits;
                for (const auto& h : hits) {
                    if (!std::isnan(h.x)) xHits.push_back(h);
                    else                   yHits.push_back(h);
                }
                clusterAlong(xHits, true);
                clusterAlong(yHits, false);
            }

            return clusters;
        }

    }
}

#endif
