#pragma once

class BarbRange {
    float hmin_, hmax_, rate_;
public:
    BarbRange(float hmin, float hmax, float rate) :
        hmin_(hmin), hmax_(hmax), rate_(rate) {}

    constexpr float hmin() const noexcept { return hmin_; }
    constexpr float hmax() const noexcept { return hmax_; }
    constexpr float rate() const noexcept { return rate_; }
    BarbRange& operator*=(float f) {
        rate_*=f;
        return (*this);
    }
    BarbRange& operator/=(float f) {
        rate_/=f;
        return (*this);
    }
};

