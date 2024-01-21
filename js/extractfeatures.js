const sampleRate = 22050;
const frameSize = 2048;
const hopLength = 512;

const midiPitches = 128;
const pitchClasses = 12;
const pitchRef = 69;
const freqRef = 440;
const windowCache = {};

let ctx = new (window.AudioContext || window.webkitAudioContext)({sampleRate: `${sampleRate}`});

// Window Function (Hann Window for example)
function applyWindowFunction(signal, windowSize) {
    if (!windowCache[windowSize]) {
        throw new Error(`Window size ${windowSize} not precomputed.`);
    }

    return signal.map((value, index) => value * windowCache[windowSize][index]);
}

// Function to zero-pad the signal to the desired length
function padSignal(signal, length) {
    if (signal.length >= length) {
        return signal.slice(0, length);
    } else {
        return [...signal, ...new Array(length - signal.length).fill(0)];
    }
}

function complexExp(theta) {
    return {
        real: Math.cos(theta),
        imag: Math.sin(theta) // Note the positive sign for the imaginary part
    };
}

// Utility functions for complex number operations
function complexAdd(a, b) {
    return { real: a.real + b.real, imag: a.imag + b.imag };
}

function complexSubtract(a, b) {
    return { real: a.real - b.real, imag: a.imag - b.imag };
}

function complexMultiply(a, b) {
    if (!a || !b) {
        console.error("Complex multiplication error: one of the operands is undefined", { a, b });
        return { real: 0, imag: 0 }; // Return a default value to avoid further errors
    }

    return {
        real: a.real * b.real - a.imag * b.imag,
        imag: a.real * b.imag + a.imag * b.real
    };
}

// Extract magnitude from FFT result (evtl. ersetzen)
function extractMagnitude(fftResult) {
    return fftResult.map(c => Math.hypot(c.real, c.imag));
}

// Function to find the nearest power of two
function nearestPowerOfTwo(num) {
    return Math.pow(2, Math.ceil(Math.log(num) / Math.log(2)));
}

// Function to zero-pad the signal to the nearest power of two
function padSignalToPowerOfTwo(signal) {
    const targetLength = nearestPowerOfTwo(signal.length);
    if (signal.length < targetLength) {
        return [...signal, ...new Array(targetLength - signal.length).fill(0)];
    }
    return signal;
}

// Precompute Hann Windows specifically for CQT
function precomputeHannWindowsForCQT(minFreq, maxFreq, binsPerOctave, sampleRate) {
    const frequencies = calculateCQTFrequencies(minFreq, maxFreq, binsPerOctave);
    frequencies.forEach(freq => {
        const windowSize = nearestPowerOfTwo(Math.floor(sampleRate / freq));
        if (!windowCache[windowSize]) {
            windowCache[windowSize] = new Array(windowSize);
            for (let i = 0; i < windowSize; i++) {
                windowCache[windowSize][i] = 0.5 * (1 - Math.cos((2 * Math.PI * i) / (windowSize - 1)));
            }
        }
    });
}

// Audio feature extraction pipeline
let analyseAudio = async (URL) => {
    ctx.resume();

    // Extract audio signal
    let signal = await x(URL);
    signal = padSignalToPowerOfTwo(signal); // Ensure signal length is a power of two

    console.log("signal:");
    console.log(signal);

    if (signal.some(isNaN)) {
        console.error("NaN found in input signal");
    }

    // Perform Constant Q Transform on the signal
    const minFreq = 20;       // Minimum frequency
    const maxFreq = 20000;    // Maximum frequency
    const binsPerOctave = 12; // Number of bins per octave

    // Precompute windows for CQT
    precomputeHannWindowsForCQT(minFreq, maxFreq, binsPerOctave, sampleRate);

    let cqtResult = constantQTransform(signal, sampleRate, minFreq, maxFreq, binsPerOctave);

    // Compute chroma features directly from CQT result
    let C = await computeChromaFromCQT(cqtResult, pitchClasses);
    // console.log("C:");
    // console.log(C);

    // Normalize chroma features
    let normC = await normalizeFeature(C, 0.001);
    // console.log("normC:");
    // console.log(normC);

    // Save some resources
    // signal = null;
    // cqtResult = null;
    // C = null;
    // ctx.suspend();

    return normC;

}

// Constant Q Transform Implementation in JavaScript

// Helper function to calculate the frequencies for each bin
function calculateCQTFrequencies(minFreq, maxFreq, binsPerOctave) {
    const factor = Math.pow(2, 1 / binsPerOctave);
    let numBins = Math.ceil(binsPerOctave * Math.log2(maxFreq / minFreq));
    let frequencies = new Array(numBins);

    for (let i = 0; i < numBins; i++) {
        frequencies[i] = minFreq * Math.pow(factor, i);
    }

    return frequencies;
}


// CQT main function
function constantQTransform(signal, sampleRate, minFreq, maxFreq, binsPerOctave) {
    const frequencies = calculateCQTFrequencies(minFreq, maxFreq, binsPerOctave);
    let cqtResult = [];

    frequencies.forEach(freq => {
        const windowSize = nearestPowerOfTwo(Math.floor(sampleRate / freq));
        const windowedSignal = applyWindowFunction(signal, windowSize);

        // Ensure the windowed signal length matches the FFT size
        const paddedWindowedSignal = padSignal(windowedSignal, windowSize);
        
        //let complexWindowedSignal = paddedWindowedSignal.map(value => ({ real: value, imag: 0 }));
        //const fftResult = cooleyTukeyFFT(complexWindowedSignal);
        
        // Usage
        let complexSignal = convertToComplexArray(signal); // Convert Float32Array to an array of complex numbers
        const fftResult = cooleyTukeyFFT(complexSignal);
        const magnitude = extractMagnitude(fftResult);

        cqtResult.push(magnitude);
    });

    console.log("cqtResult:");
    console.log(cqtResult);

    return cqtResult;
}

function cooleyTukeyFFT(signal) {
    const N = signal.length;
    if (N <= 1) return signal;

    const bitReversedSignal = bitReverseCopy(signal);

    for (let s = 1; Math.pow(2, s) <= N; s++) {
        const m = Math.pow(2, s);

        for (let k = 0; k < N; k += m) {
            let w = { real: 1, imag: 0 };

            for (let j = 0; j < m / 2; j++) {
                const theta = -2 * Math.PI * j / m;
                const wm = complexExp(theta);

                const evenIndex = k + j;
                const oddIndex = k + j + m / 2;

                // // Debugging information
                // if (oddIndex >= N) {
                //     console.error("Index out of range", { evenIndex, oddIndex, N });
                //     continue;
                // }

                const t = complexMultiply(w, bitReversedSignal[oddIndex]);
                const u = bitReversedSignal[evenIndex];

                bitReversedSignal[evenIndex] = complexAdd(u, t);
                bitReversedSignal[oddIndex] = complexSubtract(u, t);

                w = complexMultiply(w, wm);
            }
        }
    }

    return bitReversedSignal;
}

function bitReverseCopy(signal) {
    const N = signal.length;
    const bitReversedSignal = new Array(N);

    for (let i = 0; i < N; ++i) {
        const bitReversedIndex = bitReverse(i, Math.log2(N));
        bitReversedSignal[bitReversedIndex] = {...signal[i]};
    }

    return bitReversedSignal;
}

function bitReverse(value, bits) {
    let reversed = 0;
    for (let i = 0; i < bits; i++) {
        reversed = (reversed << 1) | (value & 1);
        value >>= 1;
    }
    return reversed;
}

// Compute magnitude power spectrogram, return Promise
let Y = async x => {
    console.log("Step: Magnitude spectrogram");
    let result = tf.tidy(() => tf.square(tf.abs(tf.transpose(tf.signal.stft(tf.tensor1d(x), frameSize, hopLength, frameSize, tf.signal.hannWindow)))));
    console.log(result.array());
    return result.array();
}

// Compute log-frequency spectrogram
// Y = linear frequency magnitude spectrogram
let computeSpecLogFreq = async (Y, sampleRate, frameSize, midiPitches) => {
    console.log("Step: Log-frequency spectrogram");

    // Since the number of columns isn't fixed, declaring the rows is enough (else, columns = frame count) 
    let Y_LF = new Array(midiPitches).fill(0);

    // Extract the corresponding k-th bins with k ∈ P(p) from the STFT Y
    // and sum them up
    for (let p = 0; p < midiPitches; p++) {
        let k = poolPitch(p, sampleRate, frameSize, pitchRef, freqRef);
        let lowerK = k[0];
        let upperK = k.at(-1)+1;
        let slicedY = await Y.slice(lowerK, upperK);
        Y_LF[p] = summedByCol(slicedY);
    }
    console.log(Y_LF);
    return Y_LF;
}

// Compute the chroma
// Y_LF = log-freq spectrogram
// pitchClasses = number of pitch classes, e.g. {C, C#, D, D#, E, F, F#, G, G#, A, A#, B}
let computeChroma = async (Y_LF, pitchClasses) => {
    console.log("Step: Chroma");

    let C = [];

    // Sum up all rows from the log-frequency spectrogram into their
    // respective pitch classes
    for (let c = 0; c < pitchClasses; c++) {
        C[c] = summedByCol(Y_LF.filter((e, index) => index % pitchClasses == c));
    }

    // Normalize by the max coefficient
    return await scaled(C, 1/await getMaxFromNDimArr(C));
};

// Normalize features
let normalizeFeature = async (X, threshhold) => {
    console.log("Step: Normalize");

    const K = pitchClasses;   
    const N = X[0].length;    // Frames
    let normX = arrayFilled(K, N);
    let v = 1 / Math.sqrt(K);

    for (let n = 0; n < N; n++) {

        let sumOfSquares = 0;
        for(let k = 0; k < K; k++) {
            sumOfSquares += X[k][n] ** 2;
        }
        let s = Math.sqrt(sumOfSquares);
        if (s > threshhold) {
            for(let k = 0; k < K; k++) {
                normX[k][n] = X[k][n] / s;
            }
        } else {
            for(let k = 0; k < K; k++) {
                normX[k][n] = v;
            }
        }
    }
    return normX;
}

// Compute the chroma from CQT result
// cqtResult = array of arrays, each inner array is the magnitude spectrum for a specific frequency bin
// pitchClasses = number of pitch classes (typically 12 for Western music)
let computeChromaFromCQT = async (cqtResult, pitchClasses) => {
    // Number of frames (time slices) in the CQT result
    const numFrames = cqtResult[0].length;

    // Initialize chroma matrix
    let C = arrayFilled(pitchClasses, numFrames);

    // Iterate over each frequency bin in the CQT result
    for (let bin = 0; bin < cqtResult.length; bin++) {
        // Map the CQT bin to a pitch class
        let pitchClass = bin % pitchClasses;

        // Add the magnitude spectrum of this bin to the corresponding pitch class in chroma matrix
        for (let frame = 0; frame < numFrames; frame++) {
            C[pitchClass][frame] += cqtResult[bin][frame];
        }
    }

    // Normalize by the max coefficient
    return await scaled(C, 1 / await getMaxFromNDimArr(C));
};

// Read the MIDI file
async function parseMIDI (file) {
    console.log("Step: Parse MIDI");

    let currentMidi = null;
    return new Promise((resolve, reject) => {
        const reader = new FileReader();
        reader.onload = (e) => {
        const midi = new Midi(e.target.result);
        currentMidi = midi;
        resolve(currentMidi);
    };
    reader.readAsArrayBuffer(file);
    })
}

// Compute chroma from MIDI, all in one
let analyseMIDI = async (file) => {
    console.log("Step: Analyse MIDI");

    let midi = await parseMIDI(file);

    // At ratio = 1000, there would be a column for each ms, which makes the computation needlessly slow. 
    // Thus, this reduces the temporal resolution
    let ratio = 1000 / 20;

    // This defines the time grid resolution of the features
    let C = arrayFilled(pitchClasses, Math.round(midi.duration * ratio));

    // Get duration of the entire MIDI file
    // console.log(midi);
    console.log(`--MIDI duration: ${convertHMS(midi.duration)} minutes`);

    let tracks = midi.tracks;
    let totalTracksNum = midi.tracks.length;

    for(let i = 0; i < totalTracksNum; i++) {
        let trackNotesNum = tracks[i].notes.length;

        for(let j = 0; j < trackNotesNum; j++) {
            let currentTrackNote = tracks[i].notes[j];

            // Readability
            let pitch = currentTrackNote.midi % pitchClasses;
            let noteOnTime = Math.round(currentTrackNote.time * ratio);
            let durationNote = Math.round(currentTrackNote.duration * ratio);

            // Add the notes in the feature matrix by their velocities (over all tracks)
            let noteOffTime = noteOnTime+durationNote;
            for(let k = noteOnTime; k < noteOffTime; k++) {
                C[pitch][k] = C[pitch][k] + currentTrackNote.velocity;
            }
            
            // Get info from each note
            // console.log(`
            //     MIDI: ${note.midi},
            //     Time ms: ${Math.round(note.time * 1000)},
            //     Duration ms: ${Math.round(note.duration * 1000)},
            //     Pitch: ${note.pitch},
            //     Velocity: ${note.velocity}
            //     `);
        }
    }

    // Save memory
    midi = null;

    // Return normalized chroma
    return await normalizeFeature(C, 0.001);
}