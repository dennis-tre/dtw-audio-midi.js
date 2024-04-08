function addElementWise(a, b) {
    return a.map((row, i) => row.map((val, j) => val + b[i][j]));
}

function matrixMultiply(a, b) {
    return a.map((row, i) =>
        b[0].map((_, j) =>
            row.reduce((sum, elm, k) => sum + elm * b[k][j], 0)
        )
    );
}

function getColumn(arr, columnIndex) {
    return arr.map(row => row[columnIndex]);
}

function normalizeFeatures(features, normOrder = 2, threshold = 0.001) {
    return features.map(feature => {
        const norm = feature.reduce((acc, val) => acc + Math.pow(val, normOrder), 0) ** (1 / normOrder);
        return norm > threshold ? feature.map(val => val / norm) : new Array(feature.length).fill(1 / Math.sqrt(feature.length));
    });
}

function smoothFeatures(features, windowSize = 5) {
    return features.map((feature, index, arr) => {
        let start = Math.max(0, index - Math.floor(windowSize / 2));
        let end = Math.min(arr.length, index + Math.floor(windowSize / 2) + 1);
        let sum = 0;
        for (let i = start; i < end; i++) {
            sum += arr[i];
        }
        return sum / (end - start);
    });
}

function downsampleFeatures(features, factor = 2) {
    return features.filter((_, index) => index % factor === 0);
}

function computeCostMatrix(features1, features2) {
    let costMatrix = [];
    for (let i = 0; i < features1.length; i++) {
        costMatrix[i] = [];
        for (let j = 0; j < features2.length; j++) {
            let sum = 0;
            for (let k = 0; k < features1[i].length; k++) {
                sum += Math.pow(features1[i][k] - features2[j][k], 2);
            }
            costMatrix[i][j] = Math.sqrt(sum);
        }
    }
    return costMatrix;
}

function dynamicTimeWarping(costMatrix, stepPattern) {
    let dpMatrix = Array(costMatrix.length).fill().map(() => Array(costMatrix[0].length).fill(Infinity));
    dpMatrix[0][0] = costMatrix[0][0]; // Starting point

    for (let i = 1; i < costMatrix.length; i++) {
        for (let j = 1; j < costMatrix[i].length; j++) {
            let minPreviousCost = Math.min(dpMatrix[i-1][j], dpMatrix[i][j-1], dpMatrix[i-1][j-1]);
            dpMatrix[i][j] = costMatrix[i][j] + minPreviousCost;
        }
    }

    // Backtrack to find the optimal path
    let i = costMatrix.length - 1;
    let j = costMatrix[0].length - 1;
    let path = [[i, j]];

    while (i > 0 && j > 0) {
        let steps = [[i - 1, j], [i, j - 1], [i - 1, j - 1]];
        let [nextI, nextJ] = steps.reduce((acc, step) => dpMatrix[step[0]][step[1]] < dpMatrix[acc[0]][acc[1]] ? step : acc);
        path.push([nextI, nextJ]);
        i = nextI;
        j = nextJ;
    }

    return path.reverse();
}

