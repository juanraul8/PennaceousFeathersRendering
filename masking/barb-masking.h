#pragma once
#include <list>
#include "barb-range.h"

class BarbMasking {
    float transmittance_rate_, left_rate_, right_rate_;
    std::list<BarbRange> barb_ranges_;
public:
    BarbMasking(float transmittance_rate, float left_rate, float right_rate, 
        const std::list<BarbRange>& barb_ranges) :
        transmittance_rate_(transmittance_rate), left_rate_(left_rate),
        right_rate_(right_rate),
        barb_ranges_(barb_ranges) {}
    
    constexpr float transmittance_rate() const noexcept { return transmittance_rate_; }
    constexpr float left_rate() const noexcept { return left_rate_; }
    constexpr float right_rate() const noexcept { return right_rate_; }
    const std::list<BarbRange>& barb_ranges() const noexcept { return barb_ranges_; }
};