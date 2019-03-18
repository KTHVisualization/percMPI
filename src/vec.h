#pragma once
#include <cassert>
#include <initializer_list>
#include <algorithm>
#include <array>
#include <ostream>

using ind = signed long long;

union vec3i {
    ind data[3];
    struct {
        ind x, y, z;
    };

    vec3i(ind v = 0) : x(v), y(v), z(v) {}

    vec3i(ind x, ind y, ind z) : x(x), y(y), z(z) {}

    vec3i(std::initializer_list<ind> args) {
        assert(args.size() == 3 && "Expected 3 elements.");

        auto it = args.begin();
        for (ind i = 0; i < 3; ++i) data[i] = *it++;
    }

    ind operator[](ind idx) const {
        assert(idx >= 0 && idx < 3 && "Idx out of bounds.");
        return data[idx];
    }

    ind& operator[](ind idx) {
        assert(idx >= 0 && idx < 3 && "Idx out of bounds.");
        return data[idx];
    }

    ind prod() const { return x * y * z; }

    vec3i operator+(ind s) const { return vec3i(x + s, y + s, z + s); }

    vec3i operator+(const vec3i& v) const { return vec3i(x + v.x, y + v.y, z + v.z); }
    vec3i operator-(const vec3i& v) const { return vec3i(x - v.x, y - v.y, z - v.z); }

    std::array<float, 3> operator/(const vec3i& v) const {
        return {(float)x / v.x, (float)y / v.y, (float)z / v.z};
    }

    vec3i operator*(const vec3i& v) const { return vec3i(x * v.x, y * v.y, z * v.z); }

    bool operator==(const vec3i& v) const { return x == v.x && y == v.y && z == v.z; }
    bool operator!=(const vec3i& v) const { return !(*this == v); }

    bool liesWithin(const vec3i& lower, const vec3i& upper) const {
        return x >= lower.x && x < upper.x && y >= lower.y && y < upper.y && z >= lower.z &&
               z < upper.z;
    }

    vec3i static min(const vec3i& u, const vec3i& v) {
        return vec3i(std::min(u.x, v.x), std::min(u.y, v.y), std::min(u.z, v.z));
    };
    vec3i static max(const vec3i& u, const vec3i& v) {
        return vec3i(std::max(u.x, v.x), std::max(u.y, v.y), std::max(u.z, v.z));
    };

    ind toIndexOfTotal(const vec3i& size) const { return x + y * size.x + z * size.x * size.y; }

    static vec3i fromIndexOfTotal(ind idx, const vec3i& size) {
        return vec3i(idx % (size.x), (idx / size.x) % size.y, idx / (size.x * size.y));
    }

    friend std::ostream& operator<<(std::ostream& out, const vec3i& v) {
        out << '(' << v.x << ", " << v.y << ", " << v.z << ')';
        return out;
    }
};