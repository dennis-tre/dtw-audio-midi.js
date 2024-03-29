// Calculate the DTW path for 2D time series X, Y
let dtw = async (X, Y) => {
    let n = X[0].length;
    let m = Y[0].length;
    
    // Local weights for the accumulated cost matrix D
    let w = {d: 2, h: 1.5, v: 1.5};
    console.log ("Local weights w: ", w)
      
    let C = await costMatrix(math.transpose(X), math.transpose(Y));
    let D = await computeAccumulatedCostMatrix(C, w);
    let P = await computeOptimalWarpingPath(D);

    console.log("-DTW distance DTW(X,Y) = ", D.at(-1).at(-1));

    C, R = null;

    return [P, D];
}

// Distance matrix between any 2 given 2D vectors that have the same number of rows
// Metric: Cosine distance
// Cost matrix C
let costMatrix = async (N, M) => {
    console.log("Step: Cost matrix");
    let n = N.length;
    let m = M.length;

    let cosD = arrayFilled(n, m, Infinity);

    for(let i = 0; i < n; i++) {
        for(let k = 0; k < m; k++) {
            let ni = N[i];
            let mk = M[k];
            cosD[i][k] = 1 - (math.dot(ni, mk) / (math.norm(ni, 2) * math.norm(mk, 2)));
        }
    }
    return cosD;
}


// Accumulated cost matrix D
let computeAccumulatedCostMatrix = async (C, w) => {
    console.log("Step: Accumulated cost matrix");
    let n = C.length;
    let m = C[0].length;

    // Initialize the first row and column
    for(let i = 1; i < n; i++) {
        C[i][0] = C[i-1][0] + (C[i][0] * w['h']);
    }   
    for(let j = 1; j < m; j++) {
        C[0][j] = C[0][j-1] + (C[0][j] * w['v']);
    }

    // Fill the rest of the matrix
    for(let i = 1; i < n; i++) {
        for(let k = 1; k < m; k++) {
            let C_ik = C[i][k];
            C[i][k] = Math.min( C[i-1][k] + (w['h'] * C_ik), 
                                C[i][k-1] + (w['v'] * C_ik), 
                                C[i-1][k-1] + (w['d'] * C_ik));
        }            
    }
    return C;
}


// Optimal warping path P from an accumulated cost matrix D
let computeOptimalWarpingPath = async (D) => {
    console.log("Step: Optimal warping path");
    const N = D.length;
    const M = D[0].length;

    let n = N-1;
    let m = M-1;
    let P = [[n, m]];
    let cell;

    while(n > 0 || m > 0) {
        if (n == 0) {
            cell = [0, m-1];
        } else if (m == 0) {
            cell = [n-1, 0];
        } else {
            let d_n1m1 = D[n-1][m-1];
            let d_n1m = D[n-1][m];
            let val = Math.min( d_n1m1, 
                                d_n1m, 
                                D[n][m-1] );
            if (val == d_n1m1) {
                cell = [n-1, m-1];
            } else if (val == d_n1m) {
                cell = [n-1, m];
            } else {
                cell = [n, m-1];
            }
        }
        P.push(cell);
        n = cell[0];
        m = cell[1];
    }

    P.reverse();
    return P;
}