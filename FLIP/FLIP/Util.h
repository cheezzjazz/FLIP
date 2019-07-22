#pragma once
#ifndef UTIL_H
#define UTIL_H

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

template<class T>
inline T min(const T &a1, const T &a2/*, const T &a3*/)
{
	if (a1 < a2) return a1; else return a2;
	//if (a1 < a2) return min(a1, a3); else return min(a2, a3);
}

template<class T>
inline T max(const T &a1, const T &a2/*, const T &a3*/)
{
	if (a1 < a2) return a2; else return a1;
	//if (a1 > a2) return max(a1, a3); else return max(a2, a3);
}

template<class T>
inline T sqr(const T &x)
{
	return x*x;
}

template<class T>
inline T clamp(T a, T lower, T upper)
{
	if (a < lower) return lower;
	else if (a > upper) return upper;
	else return a;
}
#endif // !UTIL_H
