// Functions that are used inside others 

// <Some overview statistics
console.log(`Sampling rate: ${sampleRate}`);
console.log(`Frame size: ${frameSize}`);
console.log(`Hop length: ${hopLength}`);
console.log(`STFT frequency resolution (in Hz): ${Math.round(((sampleRate/frameSize) + Number.EPSILON) * 100) / 100}`);
// [76, 64, 52, 40, 39, 38].forEach((value) => 
//     console.log(`Set P(${value}) = ${poolPitch(value, sampleRate, frameSize, pitchRef, freqRef)}`));


// Calculate frequency from MIDI pitch
let fPitch = (p, pitchRef, freqRef) => 2 ** ((p-pitchRef) / pitchClasses) * freqRef;

// Pool the pitches such that P(p):={k: fPitch(p-0.5) <= fPitch(k) < fPitch(p+0.5)}
let poolPitch = (p, sampleRate, frameSize, pitchRef, freqRef) => {
    
    // Calculate the border frequencies 
    let lowerPitch = fPitch(p - 0.5, pitchRef, freqRef);
    let upperPitch = fPitch(p + 0.5, pitchRef, freqRef);
    
    // Array with k-indeces ranging from 0 to (frameSize/2) + 1
    let k = [...Array((frameSize / 2) + 1).keys()];
    
    // Array with center frequency for each k
    let kFreq = k.map(x => (x * sampleRate) / frameSize);
    
    // Storing P(p) in a variable and return the corresponding k indeces
    let mask = kFreq.map(x => ((lowerPitch <= x) && (x < upperPitch)));
    return k.filter(e => mask[e]);
};

// Logarithmic compression => log(1 + gamma * v)
let logCompression = (v, gamma) => tf.log(tf.add(tf.mul(tf.tensor(gamma), v), tf.tensor(1)));

// Define a global constraint region, to which calculations are restricted
// To speed up computation. This particular configuration is called Sakoe-Chiba band
let computeGlobalConstraintRegion = async (lengthN, lengthM, windowSize) => {
    let lowerBound;
    let upperBound;
    let bw = windowSize < 1 ? ((windowSize*lengthM)/2)+1 : windowSize;
    let R = arrayFilled(2, lengthN);

    for (let i = 0; i < lengthN; i++) {
        let diagY = i/(lengthN-1)*lengthM;

        lowerBound = Math.round(diagY-bw);
        upperBound = Math.round(diagY+bw);

        if ( lowerBound < 0) lowerBound = 0;
        if ( upperBound <= 0) upperBound = bw;
        if ( lowerBound >= lengthM) lowerBound = lengthM-bw;
        if ( upperBound > lengthM) upperBound = lengthM;

        R[0][i] = lowerBound;
        R[1][i] = upperBound;
    }
    return R;
}

// Create array filled with zeros
let arrayFilled = (rows, cols, val) => {
    val = (typeof val !== 'undefined') ? val: 0;
    let filledArr = Array(rows);
    for(let i = 0; i < rows; i++) {
        filledArr[i] = Array(cols).fill(val);
    }
    return filledArr;
}

// Merge all rows from a 2D array to one summed row
let summedByCol = (arr) => arr.reduce((r, a) => {
    let al = a.length;
    for(let i = 0; i < al; i++) {
        r[i] = (r[i] || 0) + a[i];
    }
    return r;
}, []);

// Max value in 2D matrix
let getMaxFromNDimArr = async (arr) => {
    let arrL = arr.length;
    let arrPulledApart = [];
    for (let i = 0; i < arrL; i++) {
        arrPulledApart = arrPulledApart.concat(arr[i]);
    }
    let max = arrPulledApart.reduce((pv, cv) => Math.max(pv, cv), -Infinity);
    return max;
} 
// Multiply 2D array by a scalar
let scaled = async (arr, scalar) => {
    for (let r = 0; r < arr.length; r++) {
        arr[r] = arr[r].map(x => x * scalar);
    }
    return arr;
}

// Square 2D array
let square = (arr) => {
    let squaredArr = [];
    let dim0Length = arr.length;
    let dim1Length = arr[0].length;
    for (let i = 0; i < dim0Length; i++) {
        let row = [];
        for (let j = 0; j < dim1Length; j++) {
            row.push(arr[i][j] ** 2);
        }
        squaredArr.push(row);
    }
  return squaredArr;
}

let chromaHeatMap = (data, title) => {
    return {
            plot: [{
                y: ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B'],
                z: data,
                name: title,
                type: 'heatmap',
                colorscale: 'Greys',
                marker: {
                    color: '#C8A2C8',
                    line: {
                        width: 2.5
                    }
                }
            }],
            layout: {
                autosize: true,
                title: {
                    font: {
                        size: 24,
                        color: '#e8e6e3',
                    },
                    text: title,
                },
                xaxis: {
                    ticks: '',
                    side: 'bottom'
                  },
                yaxis: {
                    ticks: '',
                    ticksuffix: ' ',
                },
                paper_bgcolor: '#222',
                plot_bgcolor: '#222',
                font: {
                    color: '#eee'
                }
            },
        }
} 

let costPathPlot = (costMat, path) => {
    let layout = {
        title: {
            font: {
                size: 24,
                color: '#e8e6e3',
            },
            text: '$\\text{Accumulated cost matrix }D\\text{, }\\\\\\text{Optimal warping path }P$'
        },
        xaxis: {
            title: '$N$',
            ticks: '',
            side: 'bottom',
            gridcolor: '#575656'
        },
        yaxis: {
            title: {
                text: '$M$',
                standoff: 15,
            },
            ticks: '',
            ticksuffix: ' ',
            automargin: true,
            gridcolor: '#575656',
        },
        paper_bgcolor: '#222',
        plot_bgcolor: '#222',
        font: {
            color: '#eee'
        },
        autosize: false,
        width: 800,
        height: 800,
        margin: {
            l: 50,
            r: 50,
            b: 100,
            t: 100,
            pad: 4
        },
    };
    
    let costMatrix_D = {
        z: math.transpose(costMat),
        type: 'heatmap',
        colorscale: 'Greys',
        name: 'D',
    };

    let path_P = {
        x: path.map(e => e[0]),
        y: path.map(e => e[1]),
        type: 'scatter',
        name: 'P',
        hoverlabel: {
            font: {
                color: '#eab933'
            }
        }
    };
    
    let data = [costMatrix_D, path_P]
    
    Plotly.newPlot('path', data, layout, {scrollZoom: false});

}

function delay(time) {
    return new Promise(resolve => setTimeout(resolve, time));
  }

// For MIDI length display in console
let convertHMS = (value) => {
    const sec = parseInt(value, 10); // convert value to number if it's string
    let hours   = Math.floor(sec / 3600); // get hours
    let minutes = Math.floor((sec - (hours * 3600)) / 60); // get minutes
    let seconds = sec - (hours * 3600) - (minutes * 60); //  get seconds
    // add 0 if value < 10; Example: 2 => 02
    if (hours   < 10) {hours   = "0"+hours;}
    if (minutes < 10) {minutes = "0"+minutes;}
    if (seconds < 10) {seconds = "0"+seconds;}
    return minutes+':'+seconds; // Return is MM : SS
}