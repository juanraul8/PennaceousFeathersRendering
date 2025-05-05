#pragma once

#include <list>
#include "barb-range.h"

class FeatherMasking {
    float transmittance_rate_;
    BarbRange left_range_, right_range_;
    std::list<BarbRange> barb_ranges_;
public:
    FeatherMasking(float transmittance_rate, const BarbRange& left_range, const BarbRange& right_range, 
        const std::list<BarbRange>& barb_ranges) :
        transmittance_rate_(transmittance_rate), left_range_(left_range),
        right_range_(right_range),
        barb_ranges_(barb_ranges) {}
    
    constexpr float transmittance_rate() const noexcept { return transmittance_rate_; }
    constexpr const BarbRange& left_range() const noexcept { return left_range_; }
    constexpr const BarbRange& right_range() const noexcept { return right_range_; }
    const std::list<BarbRange>& barb_ranges() const noexcept { return barb_ranges_; }
};