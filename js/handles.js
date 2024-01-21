// Select files, execute DTW and store path data as .json
window.onload = async () => {

    let selector;
    let features = [0, 0];
    let dtwData, jsonPathData;

    // Once a file button is pressed...
    const inputFile = document.querySelectorAll('input[type=file]').forEach(e => e.addEventListener('click', (event) => {
        let elem = event.target;

        // Store pressed button as numeric value
        if (elem.id == "inputaudio") {
            selector = 0;
        } else if (elem.id == "inputmidi") {
            selector = 1;
        }

        // Retrieve window URL from input file
        elem.onchange = async () => {
            let url = null;
            const fileObj = document.querySelectorAll('input[type=file]')[selector].files[0];

            if (window.createObjcectURL != undefined) {
                url = window.createOjcectURL(fileObj);
            } else if (window.URL != undefined) {
                try {
                    url = window.URL.createObjectURL(fileObj);
                } catch (e) {
                    console.log("Abort");
                    console.error(e);
                }
            } else if (window.webkitURL != undefined) {
                url = window.webkitURL.createObjectURL(fileObj);
            }

            console.time("Extraction time");
            console.log(`__________`);

            // Once the audio file is loaded, analyse
            if (elem.id == "inputaudio") {
                features[0] = await analyseAudio(url);
                
                // When MIDI file is loaded, analyse
            } else if (elem.id == "inputmidi") {
                if (fileObj.type == "audio/mid") {
                    features[1] = await analyseMIDI(fileObj);
                    
                } else {
                    try {
                        features[1] = await analyseAudio(url);
                    }
                    catch (e) {
                        console.error("Wrong file type: ", fileObj.type, "\nAudio is needed");
                        return;
                    }
                }
            }
            console.log("Loaded: ", fileObj.name)
            console.timeEnd("Extraction time");
            console.log(`__________`);

            console.log("Current features:\n", {audio: features[0], 
                                                N: features[0][0] ? features[0][0].length : 0, 
                                                midi: features[1], 
                                                M: features[1][0] ? features[1][0].length : 0});

            // Plot features
            if (!features[0].length < 1) {
                let plotAudio = chromaHeatMap(features[0],  '$\\text{Audio STFT chroma } X_{n}$');
                Plotly.newPlot('audiofeatures', plotAudio['plot'], plotAudio['layout'], {scrollZoom: false});
            }
            if (!features[1].length < 1 ) {
                let plotMIDI = chromaHeatMap(features[1], '$\\text{MIDI chroma } Y_{m}$');
                Plotly.newPlot('midifeatures', plotMIDI['plot'], plotMIDI['layout'], {scrollZoom: false});
            }

            // Wait 5ms so the MIDI chroma has time to display
            await delay(5);

            // If both chromas are present, start calculating DTW right away
            if (!features[0] == 0 && !features[1] == 0) {
                console.log(`__________`);
                console.log("All required features are present. Starting dtw...");
                console.time("Time dtw");
                dtwData = await dtw(features[0], features[1]);
                dtwData = {path_P: dtwData[0], costMatrix_D: dtwData[1]};
                console.timeEnd("Time dtw");
                console.log(`__________`);
                console.log("Finished.\nYou can now save the dtw.json!\nData:\n", dtwData);
                
                // Plot cost matrix and path
                costPathPlot(dtwData['costMatrix_D'], dtwData['path_P'])
                jsonPathData = "text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(dtwData['path_P']));
            }
        }
    }));

    // Save and restart button actions
    const btn = Array.from(document.getElementsByClassName('btn')).forEach(e => e.addEventListener('click', async (e) => {
        let elem = e.target;

        // Empty all data
        if (elem.id == "btnReset") {
            features = [0,0];
            dtwData = null;
            jsonPathData = null;
            document.getElementById("inputaudio").value = "";
            document.getElementById("inputmidi").value = "";
            Plotly.purge("audiofeatures");
            Plotly.purge("midifeatures");
            Plotly.purge("path");
            console.clear();
            
        // Download path_P.json
        } else if (elem.id == "btnDTW") {
            if (!dtwData) { 
                console.log("Not available") 
            } else {
                console.log("Saving...");
                exportJson(jsonPathData)
            }
        }
    }));
}

// Export DTW data
function exportJson(jsonPathData) {
    keyFileName = "dtw"

    let temporaryDownloadLink = document.createElement("a");
    temporaryDownloadLink.style.display = 'none';
    document.body.appendChild( temporaryDownloadLink );

    temporaryDownloadLink.setAttribute( 'href', 'data:'+jsonPathData );
    temporaryDownloadLink.setAttribute( 'download', 'dtw.json');
    temporaryDownloadLink.click();
}