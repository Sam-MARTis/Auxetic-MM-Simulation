#pragma once
#include <cstdint>


using FlagStorageType = uint32_t;
enum class StateUpdateFlags: FlagStorageType
{
    NONE = 0,
    PARTICLE_RENDER_RADIUS = FlagStorageType{1}<<0,
    PARTICLE_RENDER_COLOUR = FlagStorageType{1}<<1,
    EDGE_RENDER_THICKNESS = FlagStorageType{1}<<2,
    EDGE_RENDER_COLOUR = FlagStorageType{1}<<3,
    RENDER_SCALE = FlagStorageType{1}<<4,
    RENDER_ORIGIN = FlagStorageType{1}<<5,
    NODE_POSITIONS = FlagStorageType{1}<<6,
    SHOULD_REBUILD_VERTICES = FlagStorageType{1}<<7,

    TIME_STEP = FlagStorageType{1}<<15,
    RELAXATION_FACTOR = FlagStorageType{1}<<16,
};

// enum class SimUpdateFlags: FlagStorageType
// {
//     NONE = 0,
//
// };

constexpr StateUpdateFlags operator|(StateUpdateFlags lhs, StateUpdateFlags rhs)
{
    return static_cast<StateUpdateFlags>(
        static_cast<FlagStorageType>(lhs) | static_cast<FlagStorageType>(rhs)
    );
}
constexpr StateUpdateFlags operator~(StateUpdateFlags flags)
{
    return static_cast<StateUpdateFlags>(
        ~static_cast<FlagStorageType>(flags)
    );
}
constexpr StateUpdateFlags operator&(StateUpdateFlags lhs, StateUpdateFlags rhs)
{
    return static_cast<StateUpdateFlags>(
        static_cast<FlagStorageType>(lhs) & static_cast<FlagStorageType>(rhs)
    );
}

constexpr StateUpdateFlags& operator&=(StateUpdateFlags& lhs, StateUpdateFlags rhs)
{
    lhs = static_cast<StateUpdateFlags>(
        static_cast<FlagStorageType>(lhs) &
        static_cast<FlagStorageType>(rhs)
    );
    return lhs;
}

constexpr StateUpdateFlags& operator|=(StateUpdateFlags& lhs, StateUpdateFlags rhs)
{
    lhs = lhs | rhs;
    return lhs;
}

inline bool any(StateUpdateFlags flags)
{
    return flags != StateUpdateFlags::NONE;
}



constexpr StateUpdateFlags RebuildRenderMask =
    StateUpdateFlags::PARTICLE_RENDER_RADIUS |
    StateUpdateFlags::PARTICLE_RENDER_COLOUR |
    StateUpdateFlags::EDGE_RENDER_THICKNESS |
    StateUpdateFlags::EDGE_RENDER_COLOUR |
    StateUpdateFlags::RENDER_SCALE |
    StateUpdateFlags::RENDER_ORIGIN |
    StateUpdateFlags::NODE_POSITIONS;

