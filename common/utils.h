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
	const vector<ui> &data;
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

	ui &operator[](ui ind)
	{
		return lookup[ind];
	}
};

class RandList
{
private:
	ui *vlist;
	ui *vpos;
	ui vnum;
	ui cap;

public:
	RandList()
	{
		vlist = vpos = nullptr;
		// vnum = cap = 0;
	};
	RandList(int _cap)
	{
		cap = _cap;
		vlist = new ui[cap];
		vpos = new ui[cap];
		vnum = 0;
		for (ui i = 0; i < cap; i++)
		{
			vpos[i] = cap;
		}
	}
	void add(int vid)
	{
		assert(vpos[vid] == cap);
		vlist[vnum] = vid;
		vpos[vid] = vnum;
		vnum++;
	};
	void remove(int vid)
	{
		assert(vpos[vid] < vnum);
		ui last_id = vlist[vnum - 1];
		ui id_pos = vpos[vid];
		vlist[id_pos] = last_id;
		vpos[last_id] = id_pos;
		vnum--;
		vpos[vid] = cap; /*set as visited*/
	}


	void clear()
	{
		for (ui i = 0; i < vnum; i++)
			vpos[i] = cap;
		vnum = 0;
	}
	ui get(ui i)
	{
		assert(i < vnum);
		return vlist[i];
	}
	ui operator[](ui i)
	{
		assert(i < vnum);
		return vlist[i];
	}
	bool contains(int vid)
	{
		return vpos[vid] != cap;
	}
	bool empty() { return vnum == 0; }
	ui size() { return vnum; }
	ui getCap() { return cap; }
	void dispose()
	{
		if (vlist != nullptr)
		{
			delete[] vlist;
			vlist = nullptr;
		}
		if (vpos != nullptr)
		{
			delete[] vpos;
			vpos = nullptr;
		}
	}
	~RandList()
	{
		dispose();
	}
#ifdef DBGMOD
	void printList(FILE *f = stdout)
	{
		fprintf(f, "Total %d: ", vnum);
		int *tmp_lst = new int[cap];
		memcpy(tmp_lst, vlist, vnum * sizeof(int));
		std::sort(tmp_lst, tmp_lst + vnum);
		// qsort(tmp_lst, vnum, sizeof(int), cmpfunc);
		for (ui i = 0; i < vnum; i++)
		{
			fprintf(f, "%d ", tmp_lst[i]);
		}
		fprintf(f, "\n");
	};
#else
	void printList(FILE *f = stdout){};
#endif
};