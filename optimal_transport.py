from functools import reduce
from collections import deque

def show(arr):
    return '\n'.join('  | '.join(map(str, row)) for row in arr)
        
class Pair:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        
    def __str__(self):
        return f"{self.x} {'+' if self.y >= 0 else '-'} {abs(self.y)}·ε"
        
    def __neg__(self, other):
        return Pair(-self.x, -self.y)
        
    def is_pos(self):
        return self > Pair(0, 0)
        
    def __add__(self, other):
        return Pair(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return Pair(self.x - other.x, self.y - other.y)

    def __lt__(self, other):
        return (self.x, self.y) < (other.x, other.y)

    def __le__(self, other):
        return (self.x, self.y) <= (other.x, other.y)      

    def __gt__(self, other):
        return (self.x, self.y) > (other.x, other.y)

    def __ge__(self, other):
        return (self.x, self.y) > (other.x, other.y)  

    def __mul__(self, k):
        assert isinstance(k, int), "Multiplication undefined"
        return Pair(k * self.x, k * self.y)

    def __rmul__(self, k):
        assert isinstance(k, int), "Multiplication undefined"
        return Pair(k * self.x, k * self.y)
    

def pair_sum(xs):
    s = Pair(0, 0)
    for x in xs:
        s += x
    return s

        
class Simplex:
    def __init__(self, a, b, c):
        self.m, self.n = len(a), len(b)
        self.costs = c
        self.a = [Pair(x, 1) for i, x in enumerate(a)]
        self.b = [Pair(x, self.m) if i == self.n - 1 else Pair(x, 0) for i, x in enumerate(b)]
              
    def make_initial_solution(self):
        inv_pos = {}
        pos, row, col = [], {i: set() for i in range(self.m)}, {j: set() for j in range(self.n)}, 
        
        
        arr = [[Pair(0, 0) for _ in range(self.n + 1)] for _ in range(self.m + 1)]
        for i in range(self.m):
            arr[i][self.n] = self.a[i]
        for j in range(self.n):
            arr[self.m][j] = self.b[j]
        i, j = 0, 0
        num_pos = 0

        while i < self.m or j < self.n: 
            sum_row = arr[i][self.n] - pair_sum([arr[i][k] for k in range(j)]) if i < self.m else Pair(0, 0)
            sum_col = arr[self.m][j] - pair_sum([arr[k][j] for k in range(i)]) if j < self.n else Pair(0, 0)                                                                         
            if not sum_row.is_pos():
                i += 1    
            elif not sum_col.is_pos():
                j += 1
            else:
                if sum_row < sum_col:
                    arr[i][j] = sum_row
                    pos.append((i, j))
                    inv_pos[(i, j)] = num_pos
                    num_pos += 1
                    row[i].add((i, j))
                    col[j].add((i, j))
                    i += 1                    
                else:
                    arr[i][j] = sum_col
                    pos.append((i, j))
                    inv_pos[(i, j)] = num_pos
                    num_pos += 1
                    row[i].add((i, j))
                    col[j].add((i, j))
                    j += 1

            arr[self.m][self.n] = pair_sum([arr[self.m][k] for k in range(self.n)])
            
        k = len(pos)
        N = self.n + self.m - 1
        if k < N:
            return arr, pos, inv_pos, row, col, True
   

        return arr, pos, inv_pos, row, col, False


    def step(self, arr, pos, inv_pos, row, col):
        x, y = pos[0]
        u, v, visited = [None for _ in range(self.m)], [None for _ in range(self.n)], set([(x, y)])
        lam, mu = {i: {} for i in range(self.m)}, {i: {} for i in range(self.n)}
        lam[x][0] = 1
        u[x] = self.costs[x][y]
        v[y] = 0
        _pos = set(pos[1:])
        k = 0
        while _pos:
            if k > 0:
                x, y = _pos.pop()
                
            queue = deque([(x, y)])
            while queue:
                i, j = queue.pop()
                for a, b in row[i]:
                    if not (a, b) in visited and (a, b) != (i, j) and v[b] is None:
                        v[b] = v[j] + self.costs[a][b] - self.costs[i][j]
                        i1, i2 = inv_pos[(i, j)], inv_pos[(a, b)]
                        d = dict(mu[j])
                        d[i2] = d.get(i2, 0) + 1
                        d[i1] = d.get(i1, 0) - 1
                        mu[b] = d
                        visited.add((a, b))
                        queue.append((a, b))
                        _pos.discard((a, b))
                        k += 1
                for a, b in col[j]:
                    if not (a, b) in visited  and (a, b) != (i, j) and u[a] is None:
                        u[a] = u[i] + self.costs[a][b] - self.costs[i][j]
                        i1, i2 = inv_pos[(i, j)], inv_pos[(a, b)]
                        d = dict(lam[i])
                        d[i2] = d.get(i2, 0) + 1
                        d[i1] = d.get(i1, 0) - 1
                        lam[a] = d
                        visited.add((a, b))
                        queue.append((a, b))
                        _pos.discard((a, b))
                        k += 1

        max_delta = 0
        for i in range(self.m):
            for j in range(self.n):
                delta = u[i] + v[j] - self.costs[i][j]
                if delta > max_delta:
                    max_delta = delta
                    i0, j0 = i, j

            return (arr, pos, inv_pos, row, col, max_delta)

        ##############
        theta = Pair(float('inf'), float('inf'))
        theta_idx = set([(-1, (i0, j0))])
        x_rem, y_rem = None, None
        for k in range(self.n + self.m - 1):
            nu = lam[i0].get(k, 0) + mu[j0].get(k, 0)
            x, y = pos[k]
            if nu == 1:
                if arr[x][y] < theta:
                    theta = arr[x][y]
                    x_rem, y_rem = x, y
                theta_idx.add((1, (x, y)))
            elif nu == -1:
                theta_idx.add((-1, (x, y)))
        if x_rem is None:
            return (arr, pos, inv_pos, row, col, 0)
        for v, (x, y) in theta_idx:
            arr[x][y] -= v * theta
        
        idx = inv_pos[(x_rem, y_rem)]
        pos[idx] = (i0, j0)
        del inv_pos[(x_rem, y_rem)]
        inv_pos[(i0, j0)] = idx
        row[x_rem].remove((x_rem, y_rem))
        col[y_rem].remove((x_rem, y_rem))
        row[i0].add((i0, j0))
        col[j0].add((i0, j0))
        return arr, pos, inv_pos, row, col, max_delta 

    def solve(self):
        arr, pos, inv_pos, row, col, done = self.make_initial_solution()
        if not done:
            max_delta = float('inf')
            k = 0
            while max_delta > 0:
                arr, pos, inv_pos, row, col, max_delta = self.step(arr, pos, inv_pos, row, col)
                k += 1
        ans = 0
        for i in range(self.m):
            for j in range(self.n):
                ans += self.costs[i][j] * arr[i][j].x
        return ans
    
def minimum_transportation_price(a, b, c):
    simplex= Simplex(a, b, c)
    return simplex.solve()
