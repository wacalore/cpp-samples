#include <iostream>

using namespace std;


char floor3[3][4] =
{
"1..",
"oo.",
"...",
};

char floor2[3][4] =
{
"ooo",
"..o",
".oo",
};

char floor1[3][4] =
{
"ooo",
"o..",
"o.2",
};

template <size_t n, size_t size_y> int findPrincess(int x, int y, char (&floor)[n][size_y])
{
    int m = 1000; int p=0;
    if (floor[x][y] == '2') return 0;

    floor[x][y] = 'x';

    if (floor == floor3)
    {
        if(floor2[x][y] == '.') {m = 5 + findPrincess(x,y,floor2);}
        if(x+1 < n+1 && floor[x+1][y] == '.')
        {
            p = 5 + findPrincess(x+1,y,floor);
            if (p < m) { m = p;}
        }
        if(y+1 < size_y && floor[x][y+1] == '.')
        {
            p = 5 + findPrincess(x, y+1, floor);
            if (p < m) { m = p;}
        }
        if(x-1 >= 0 && floor[x-1][y] == '.')
        {
            p = 5 + findPrincess(x-1, y, floor);
            if (p < m) { m = p; }
        }
        if(y-1 >= 0 && floor[x][y-1] == '.')
        {
            p = 5 + findPrincess(x, y-1, floor);
            if (p < m) { m = p; }
        }
        return m;
    }

    if (floor == floor2)
    {
        if(floor1[x][y] == '.')
        {
            m = 5 + findPrincess(x,y,floor1);
        }
        if(x+1 < n+1 && floor[x+1][y] == '.')
        {
            p = 5 + findPrincess(x+1,y,floor);
            if (p < m) { m = p;}
        }
        if(y+1 < size_y && floor[x][y+1] == '.')
        {
            p = 5 + findPrincess(x, y+1, floor);
            if (p < m) { m = p;}
        }
        if(x-1 >= 0 && floor[x-1][y] == '.')
        {
            p = 5 + findPrincess(x-1, y, floor);
            if (p < m) { m = p; }
        }
        if(y-1 >= 0 && floor[x][y-1] == '.')
        {
            p = 5 + findPrincess(x, y-1, floor);
            if (p < m) { m = p; }
        }
        return m;
    }

    if(floor == floor1)
    {
        if(x+1 < n+1 && (floor[x+1][y] == '.' || floor[x+1][y] == '2'))
        {
            m = 5 + findPrincess(x+1,y,floor);
        }
        if(y+1 < size_y && (floor[x][y+1] == '.' || floor[x][y+1] == '2'))
        {
            p = 5 + findPrincess(x, y+1, floor);
            if (p < m) { m = p;}
        }
        if(x-1 >= 0 && (floor[x-1][y] == '.' || floor[x-1][y] == '2'))
        {
            p = 5 + findPrincess(x-1, y, floor);
            if (p < m) { m = p; }
        }
        if(y-1 >= 0 && (floor[x][y-1] == '.' || floor[x][y-1] == '2'))
        {
            p = 5 + findPrincess(x, y-1, floor);
            if (p < m) { m = p; }
        }
        return m;
    }

}

int main()
{
    cout << "Found Princess in : " << findPrincess(0, 0, floor3) << " seconds." << endl;
    return 0;
}
