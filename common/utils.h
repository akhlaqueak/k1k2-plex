#include "graph.h"
class MBitSet
{
private:
public:
	ui n;
	ui cap;
	ui *buf;

	MBitSet()
	{
		buf = nullptr;
		cap = n = 0;
	}
	MBitSet(ui _cap)
	{
		cap = _cap;
		n = (cap >> 5) + 1;
		buf = new ui[n];
		fill(buf, buf + n, 0);
		// for (ui i = 0; i < n; ++i)
		// 	buf[i] = 0;
	}
	~MBitSet()
	{
		// todo doing double free, see what causing it and fix

		if (buf != nullptr)
			// delete[] buf;
			buf = nullptr;
	}
	void reset()
	{
		fill(buf, buf + n, 0);
	}

	void setAll()
	{
		fill(buf, buf + n, 0xffffffff);
	}
	// FLIP all the bits
	void flip()
	{
		for (ui i = 0; i < n; ++i)
			buf[i] = ~buf[i];
	}
	void set(ui x)
	{
		// assert(x < cap);
		buf[x >> 5] |= (ui)1 << (x & 31);
	}

	bool test(ui x)
	{
		// cout << x << " " << n << " " << cap << endl;
		return buf[x >> 5] >> (x & 31) & 1;
	}

	bool empty()
	{
		for (ui i = 0; i < n; ++i)
			if (buf[i])
				return false;
		return true;
	}
};

class Lookup
{
	vector<ui> &lookup;
	const vector<ui>& data;
	bool binary;

public:
	Lookup(vector<ui> &_lookup, vector<ui> &_data, bool mode = false) : lookup(_lookup), data(_data)
	{
		binary = mode;
		for (ui ind = 0; ind < data.size(); ind++)
		{
			ui u = data[ind];
			if (binary)
				lookup[u] = 1;
			else
				lookup[u] = ind + 1;
		}
	}

	~Lookup()
	{
		for (const ui &u : data)
		{
			lookup[u] = 0;
		}
	}

	void erase()
	{
		for (const ui &u : data)
		{
			lookup[u] = 0;
		}
	}

	ui& operator[](ui ind)
	{
		return lookup[ind];
	}
};