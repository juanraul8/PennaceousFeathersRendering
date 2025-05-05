#pragma once

class BarbuleMasking {
    float transmittance_rate_;
    float hmin_, hmax_;
public:
    BarbuleMasking(float transmittance_rate, float hmin = -1.0f, float hmax = 1.0f) :
        transmittance_rate_(transmittance_rate), hmin_(hmin), hmax_(hmax) {}
    
    constexpr float transmittance_rate() const noexcept { return transmittance_rate_; }
    constexpr float barbule_rate() const noexcept { return 1.0f - transmittance_rate(); }
    constexpr float hmin() const noexcept { return hmin_; }
    constexpr float hmax() const noexcept { return hmax_; }
};