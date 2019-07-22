#pragma once

#include <cassert>
#include <iostream>
#include <cmath>

template<class T>
struct Vec2
{
	T v[2];

	Vec2() {};
	
	Vec2(T value_)
	{
		v[0] = v[1] = value_;
	}

	Vec2(T v0, T v1)
	{
		v[0] = v0; v[1] = v1;
	}

	T &operator[](int index)
	{
		assert(0 <= index && (unsigned int)index<2);
		return v[index];
	}

	const T &operator[](int index) const
	{
		assert(0 <= index && (unsigned int)index<2);
		return v[index];
	}

	Vec2<T> operator+=(const Vec2<T> &w)
	{
		v[0] += w.v[0]; v[1] += w.v[1]; return *this;
	}

	template<class S>
	Vec2<T> operator*=(S scalar)
	{
		v[0] *= scalar; v[1] *= scalar; return *this;
	}
};

template<class S, class T>
inline Vec2<T> operator*(const Vec2<T> &a, S scalar)
{
	return Vec2<T>(scalar*a.v[0], scalar*a.v[1]);
}

template<class S, class T>
inline Vec2<T> operator*(S scalar, const Vec2<T> &a)
{
	return Vec2<T>(scalar*a.v[0], scalar*a.v[1]);
}

template<class T> inline Vec2<T> operator+(const Vec2<T> &a, const Vec2<T> &b)
{
	return Vec2<T>(a.v[0] + b.v[0], a.v[1] + b.v[1]);
}