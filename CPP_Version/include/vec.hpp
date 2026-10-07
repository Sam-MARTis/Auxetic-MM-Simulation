#pragma once

#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <cassert>
#include <array>

template<typename T, uint N>
struct Vec
{
    std::array<T, N> data;
    Vec()
    {
        data.fill(T{});
    }
    explicit Vec(T val)
    {
        data.fill(val);
    }
    Vec(std::initializer_list<T> values)
    {
        // data.fill(T{});
        assert(N == values.size());
        std::copy(values.begin(), values.end(), data.begin());
    }
    T& operator[](uint i)
    {
        return data[i];
    }
    const T& operator[](size_t i) const
    {
        return data[i];
    }
    void operator+=(const Vec& vec)
    {
        for (uint i = 0; i < N; i++)
        {
            data[i] += vec.data[i];
        }
    }
    void operator-=(const Vec& vec)
    {
        for (uint i = 0; i < N; i++)
        {
            data[i] -= vec.data[i];
        }
    }
    Vec operator+(const Vec& vec) const
    {
        vec result(*this);
        result += vec;
        return result;
    }
    Vec operator-(const Vec& vec) const
    {
        vec result(*this);
        result -= vec;
        return result;
    }
    void operator*=(const T val)
    {
        for (uint i = 0; i < N; i++)
        {
            data[i] *= val;
        }
    }
    void operator/=(const T val)
    {
        *this *= 1/val;
    }
    Vec operator*(const T val) const
    {
        Vec result(*this);
        result *= val;
        return result;
    }
    friend Vec operator*(T scalar, const Vec& v) {
        return v * scalar;
    }

    Vec operator/(const T val) const
    {
        Vec result(*this);
        result /= val;
        return result;
    }
    [[nodiscard]] T norm() const
    {
        T sum = 0;
        for (uint i = 0; i < N; i++)
        {
            sum += data[i] * data[i];
        }
        return std::sqrt(sum);
    }
    constexpr static T dot(const Vec& a, const Vec& b)
    {
        T sum = 0;
        for (uint i = 0; i < N; i++)
        {
            sum += a.data[i] * b.data[i];
        }
        return sum;
    }
    constexpr static T dot(const Vec& a, const std::array<T, N>& b)
    {
        T sum = 0;
        for (uint i = 0; i < N; i++)
        {
            sum += a.data[i] * b[i];
        }
        return sum;
    }
};

namespace VecMath
{
    template <typename T, uint N, uint M>
    Vec<T, N> matrixMultiply(const Vec<Vec<T, M>, N>& matrix, const Vec<T, M>& vector)
    {
        Vec<T, N> result;
        for (uint i = 0; i < N; ++i)
        {
            result[i] = Vec<T, M>::dot(matrix[i], vector);
        }
        return result;
    }

    template <typename T, uint N>
    Vec<T, N> componentwiseMultiply(const Vec<T, N>& a, const Vec<T, N>& b)
    {
        Vec<T, N> result;
        for (uint i = 0; i < N; ++i)
        {
            result[i] = a[i] * b[i];
        }
        return result;
    }

    inline bool orientation(const Vec<double, 2>& p, const Vec<double, 2>& q, const Vec<double, 2>& r)
    {
        double val = (q[1]-p[1])*(r[0] - q[0]) - (q[0]-p[0])*(r[1] - q[1]);
        return val <= 0;
    }
    inline bool segmentsIntersects2D(const Vec<double, 2>& e1p1, const Vec<double, 2>& e1p2, const Vec<double, 2>& e2p1, const Vec<double, 2>& e2p2)
    {
        const bool o1 = orientation(e1p1, e1p2, e2p1);
        const bool o2 = orientation(e1p1, e1p2, e2p2);
        const bool o3 = orientation(e2p1, e2p2, e1p1);
        const bool o4 = orientation(e2p1, e2p2, e1p2);
        return o1!=o2 && o3!=o4;
    }

}



