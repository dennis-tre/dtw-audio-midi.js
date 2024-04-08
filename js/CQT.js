function kron(a, b) {
    let aRows = a.length, aCols = a[0].length;
    let bRows = b.length, bCols = b[0].length;
    let result = new Array(aRows * bRows);

    for (let i = 0; i < result.length; i++) {
        result[i] = new Array(aCols * bCols);
    }

    for (let i = 0; i < aRows; i++) {
        for (let j = 0; j < aCols; j++) {
            for (let k = 0; k < bRows; k++) {
                for (let l = 0; l < bCols; l++) {
                    result[i * bRows + k][j * bCols + l] = a[i][j] * b[k][l];
                }
            }
        }
    }

    return result;
}

function flatten(array) {
    return array.reduce((acc, val) => acc.concat(Array.isArray(val) ? flatten(val) : val), []);
}

function reshape(flatArray, shape) {
    function _reshape(subArray, subShape) {
        if (subShape.length === 1) {
            return subArray.splice(0, subShape[0]);
        }
        const len = subShape.shift();
        let newArray = [];
        for (let i = 0; i < len; i++) {
            newArray.push(_reshape(subArray, subShape.slice()));
        }
        return newArray;
    }
    return _reshape(flatArray.slice(), shape.slice());
}

function getShape(array) {
    let shape = [];
    while (Array.isArray(array)) {
        shape.push(array.length);
        array = array[0];
    }
    return shape;
}

function roll(array, shift, axis = null) {
    if (axis === null) {
        // Flatten, roll, and reshape
        const flatArray = flatten(array);
        const shape = getShape(array);
        const rolledFlatArray = roll(flatArray, shift);
        return reshape(rolledFlatArray, shape);
    } else {
        // Roll along a specific axis
        // This simple example assumes rolling along the first axis for simplicity
        let result = [];
        const len = array.length;
        for (let i = 0; i < len; i++) {
            result[(i + shift + len) % len] = array[i];
        }
        return result;
    }
}

// This implementation provides a basic framework. It handles flattening and reshaping for a generic roll operation
// but rolling along a specific axis other than the first requires additional logic based on the array's dimensions and the specified axis.

////

function rollAxis(array, shift, axis) {
    const shape = getShape(array);
    if (axis < 0 || axis >= shape.length) {
        throw new Error("Axis out of bounds");
    }

    function recursiveRoll(arr, currentAxis) {
        if (currentAxis === axis) {
            // This is the axis we want to roll along, perform the shift here
            const len = arr.length;
            const result = new Array(len);
            for (let i = 0; i < len; i++) {
                result[(i + shift + len) % len] = arr[i];
            }
            return result;
        } else {
            // Recurse into the next dimension
            return arr.map(subArr => recursiveRoll(subArr, currentAxis + 1));
        }
    }

    return recursiveRoll(array, 0);
}

// Helper functions `flatten`, `reshape`, `getShape` remain as previously defined

////

function concatenateArrays(...arrays) {
    return arrays.reduce((acc, curr) => acc.concat(curr), []);
}

function createRange(start, stop, step = 1) {
    let result = [];
    for (let i = start; i < stop; i += step) {
        result.push(i);
    }
    return result;
}

function parseSlice(sliceString) {
    let [start, stop, step] = sliceString.split(':').map(Number);
    // Default values if undefined
    start = isNaN(start) ? 0 : start;
    step = isNaN(step) ? 1 : step;

    if (isNaN(stop)) {
        throw new Error("Stop value is required in slice notation.");
    }

    return createRange(start, stop, step);
}

// Adjusting `createRange` to handle negative steps for reverse sequences
function createRange(start, stop, step = 1) {
    let result = [];
    if (step > 0) {
        for (let i = start; i < stop; i += step) {
            result.push(i);
        }
    } else {
        for (let i = start; i > stop; i += step) { // Handle negative steps
            result.push(i);
        }
    }
    return result;
}

function upsample(s, n, phase = 0) {
    // Create the sequence with 1 followed by n-1 zeros
    let sequence = [1].concat(Array(n-1).fill(0));
    
    // Use the previously defined kron function for the Kronecker product
    let kronResult = kron(s, sequence);
    
    // Use the previously defined roll function to shift by phase
    let result = roll(kronResult, phase);
    
    return result;
}

function buffer(x, n, p = 0, opt = null) {
    if (opt !== 'nodelay' && opt !== null) {
        throw new Error(`${opt} not implemented`);
    }

    let result = [];
    let i = 0;

    if (opt === 'nodelay') {
        // Directly start filling the buffer without initial zeros
        i = 0;
    } else {
        // Optionally start with `p` zeros
        // This handles cases where the initial buffer is partially filled with zeros and then x
        let zeros = p > 0 ? new Array(p).fill(0) : [];
        let firstFrame = zeros.concat(x.slice(0, n - zeros.length));
        result.push(firstFrame);
        i = firstFrame.length - zeros.length;
    }

    while (i < x.length) {
        let frame;
        if (i + n - p <= x.length) {
            // Normal frame
            frame = x.slice(i, i + n - p);
        } else {
            // Last frame, potentially shorter than n
            frame = x.slice(i, x.length);
        }

        if (p > 0 && result.length > 0) {
            // Add `p` elements from the end of the previous frame to the start of the new frame
            frame = result[result.length - 1].slice(n - p).concat(frame);
        }

        if (frame.length < n) {
            // If the frame is shorter than n, pad it with zeros at the end
            frame = frame.concat(new Array(n - frame.length).fill(0));
        }

        result.push(frame);
        i += n - p;
    }

    // Transpose the result to match the Python and MATLAB output
    if (result.length > 0) {
        let transposedResult = [];
        for (let col = 0; col < n; col++) {
            transposedResult[col] = result.map(row => row[col]);
        }
        return transposedResult;
    }

    return result;
}

