
#pragma once

#include "mesh.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef SHOULD_PROFILE
namespace profiling
{
    struct ProfileStats
    {
        std::uint64_t calls = 0;
        double totalNs = 0.0;
        double minNs = std::numeric_limits<double>::max();
        double maxNs = 0.0;

        void add(const double elapsedNs)
        {
            calls++;
            totalNs += elapsedNs;
            minNs = std::min(minNs, elapsedNs);
            maxNs = std::max(maxNs, elapsedNs);
        }
    };

    class Registry
    {
    public:
        static Registry& instance()
        {
            static Registry registry;
            return registry;
        }

        void add(const std::string& label, const double elapsedNs)
        {
            samples[label].add(elapsedNs);
        }

        void print(std::ostream& os = std::cout) const
        {
            if (samples.empty())
            {
                os << "[profiling] No samples recorded.\n";
                return;
            }

            struct Row
            {
                std::string label;
                ProfileStats stats;
            };

            std::vector<Row> rows;
            rows.reserve(samples.size());
            for (const auto& [label, stats] : samples)
            {
                rows.push_back(Row{label, stats});
            }

            std::sort(rows.begin(), rows.end(), [](const Row& lhs, const Row& rhs)
            {
                return lhs.stats.totalNs > rhs.stats.totalNs;
            });

            const auto nsToMs = [](const double ns)
            {
                return ns / 1'000'000.0;
            };

            std::size_t labelWidth = std::string("Scope").size();
            for (const auto& row : rows)
            {
                labelWidth = std::max(labelWidth, row.label.size());
            }

            const int callsWidth = 10;
            const int timeWidth = 12;

            os << "\nProfiling summary (times in ms)\n";
            os << std::left << std::setw(static_cast<int>(labelWidth)) << "Scope"
               << " | " << std::right << std::setw(callsWidth) << "Calls"
               << " | " << std::setw(timeWidth) << "Total"
               << " | " << std::setw(timeWidth) << "Average"
               << " | " << std::setw(timeWidth) << "Min"
               << " | " << std::setw(timeWidth) << "Max"
               << '\n';

            os << std::string(labelWidth, '-') << "-+-"
               << std::string(callsWidth, '-') << "-+-"
               << std::string(timeWidth, '-') << "-+-"
               << std::string(timeWidth, '-') << "-+-"
               << std::string(timeWidth, '-') << "-+-"
               << std::string(timeWidth, '-') << '\n';

            os << std::fixed << std::setprecision(3);
            for (const auto& row : rows)
            {
                const double averageNs = row.stats.calls == 0 ? 0.0 : row.stats.totalNs / static_cast<double>(row.stats.calls);
                const double minNs = row.stats.calls == 0 ? 0.0 : row.stats.minNs;
                const double maxNs = row.stats.calls == 0 ? 0.0 : row.stats.maxNs;

                os << std::left << std::setw(static_cast<int>(labelWidth)) << row.label
                   << " | " << std::right << std::setw(callsWidth) << row.stats.calls
                   << " | " << std::setw(timeWidth) << nsToMs(row.stats.totalNs)
                   << " | " << std::setw(timeWidth) << nsToMs(averageNs)
                   << " | " << std::setw(timeWidth) << nsToMs(minNs)
                   << " | " << std::setw(timeWidth) << nsToMs(maxNs)
                   << '\n';
            }
        }

    private:
        std::unordered_map<std::string, ProfileStats> samples;
    };

    class ScopeTimer
    {
    public:
        explicit ScopeTimer(const std::string& label)
            : label(label), start(Clock::now())
        {
        }

        ~ScopeTimer()
        {
            const auto elapsed = Clock::now() - start;
            const double elapsedNs = std::chrono::duration<double, std::nano>(elapsed).count();
            Registry::instance().add(label, elapsedNs);
        }

    private:
        using Clock = std::chrono::steady_clock;

        std::string label;
        Clock::time_point start;
    };

    inline void printReport(std::ostream& os = std::cout)
    {
        Registry::instance().print(os);
    }

    inline void printLiveReport(std::ostream& os = std::cout)
    {
        os << "\x1b[2J\x1b[H";
        Registry::instance().print(os);
        os.flush();
    }
}

#define PROFILE_CONCAT_INNER(x, y) x##y
#define PROFILE_CONCAT(x, y) PROFILE_CONCAT_INNER(x, y)
#define PROFILE_SCOPE(label) ::profiling::ScopeTimer PROFILE_CONCAT(_profileScope_, __COUNTER__){label}
#else
namespace profiling
{
    inline void printReport(std::ostream& os = std::cout)
    {
        (void)os;
    }
}

#define PROFILE_SCOPE(label) do { (void)sizeof(label); } while (false)
#endif

struct BoundaryConstraint
{
    Vec2D nVec;
    double c;
    unsigned int boundaryId;
};

class Solid
{
    SpMat A;
    Vec c;
    SpMat K;
    Vec g;
private:
    std::vector<BoundaryConstraint> constraints;

    // Remember to call flat vec initialization of mesh when solid made
    // void intializeNormalDirichletBCs();
    void assembleConstraintMatrix();
    void buildConstraintRHS();
    void assembleGradientAndHessian();
    void solve();
public:
    Solid(Mesh& solidMesh);
    Mesh& mesh;

    void timeStep(const double dt);

};
