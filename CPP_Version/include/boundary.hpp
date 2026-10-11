#pragma once

#include <cstdint>

#include "constants.hpp"
#include "utils.hpp"
#include "vec.hpp"

// Boundary conditions for the static energy minimiser; see DeveloperReference.md
// section 3 for how each type maps onto the solver.
enum class BoundaryConditionType : std::uint32_t
{
    FREE = 0,
    FIXED,
    FIXED_NORMAL,
    DISPLACE_NORMAL,
};

struct BoundaryCondition
{
    BoundaryConditionType type = BoundaryConditionType::FREE;

    // Signed rate along +n; positive is normally outward. DISPLACE_NORMAL only.
    double velocity = 0.0;

    // Explicit unit normal; (0,0) derives it from the tagged geometry.
    Vec2D normalOverride{0.0, 0.0};

    [[nodiscard]] bool usesNormalConstraint() const
    {
        return type == BoundaryConditionType::FIXED_NORMAL ||
               type == BoundaryConditionType::DISPLACE_NORMAL;
    }

    static BoundaryCondition free()
    {
        return BoundaryCondition{};
    }

    static BoundaryCondition fixed()
    {
        BoundaryCondition condition;
        condition.type = BoundaryConditionType::FIXED;
        return condition;
    }

    static BoundaryCondition fixedNormal(const Vec2D& normalOverride = Vec2D{0.0, 0.0})
    {
        BoundaryCondition condition;
        condition.type = BoundaryConditionType::FIXED_NORMAL;
        condition.normalOverride = normalOverride;
        return condition;
    }

    static BoundaryCondition displaceNormal(const double velocity,
                                            const Vec2D& normalOverride = Vec2D{0.0, 0.0})
    {
        BoundaryCondition condition;
        condition.type = BoundaryConditionType::DISPLACE_NORMAL;
        condition.velocity = velocity;
        condition.normalOverride = normalOverride;
        return condition;
    }
};

// A boundary chord; per-node normals are derived from these.
struct BoundarySegment
{
    Index node1 = NO_BOUNDARY;
    Index node2 = NO_BOUNDARY;
    Index boundaryId = NO_BOUNDARY;
    Vec2D outwardNormal{0.0, 0.0};
};

inline const char* boundaryConditionTypeName(const BoundaryConditionType type)
{
    switch (type)
    {
        case BoundaryConditionType::FREE: return "FREE";
        case BoundaryConditionType::FIXED: return "FIXED";
        case BoundaryConditionType::FIXED_NORMAL: return "FIXED_NORMAL";
        case BoundaryConditionType::DISPLACE_NORMAL: return "DISPLACE_NORMAL";
    }
    return "UNKNOWN";
}
