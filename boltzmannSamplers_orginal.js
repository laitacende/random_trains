const {loga, geometric, bernoulli, poisson, logaTruncated, geometricTruncated} = require("./distributions.js");
const fs = require('node:fs');
let endpointPlank = "○";
let endpointWheel = "•";
let endpointPerson = "p"

function endpointEGF(x) {
    return x; // rozmiar 1,  tworząca == 1 * x/1! - jeden element wielkości 1 i jedno możliwe wartościowanie
}

function gammaEndpoint(x, symbol) {
    return symbol;
}

function gammaCycEndpoint(x, symbol, start) {
    let k = loga(endpointEGF(x), start);
    if (k == -1) {
        return -1;
    }
    // make k copies of GammaEndpoint(endpointGF))
    cyc = []
    for (let i = 0; i < k; i++) {
        cyc.push(gammaEndpoint(endpointEGF(x), symbol));
    }
    return cyc;
}

function cycEGF(x) {
    return Math.log(1.0 / (1.0 - x));
}

function cycDEGF(x) {
    return 1.0 / (1 - x);
}

function cycDDEGF(x) {
    return 1.0 / ((1.0 - x) * (1.0 - x));
}


function gammaDisjointUnion(x, symbol, start) {
    let p = 1.0 / (1.0 + cycEGF(x)); // Uwaga: prawdopodobieństwo, że kół nie będzie jest bardzo duże...
                                // maleje wraz z wzrostem x, ale później taki x jest za duży dla zbieżności train...
    let tmp = bernoulli(p);
    if (tmp == 1) {
        return [];
    } else {
        let tmp1 = gammaCycEndpoint(x, symbol, 1);
        if (tmp1 != -1) {
            return [tmp1];
        } else {
            return -1;
        }
    }
}


// Z * Z * (1 + Cyc(Z))
function gammaPlanks(x) {
    // two endpoints
    let plank = [gammaEndpoint(x, endpointPlank),
            gammaEndpoint(x, endpointPlank)];
    let res = gammaDisjointUnion(x, endpointWheel, 1);
    if (res == -1) {
        return -1;
    }
    return {"planks": plank, "wheels": res};
}

// x^2(1+log(1 / (1-x)))
function plankEGF(x) {
    return x * x * (1 + cycEGF(x));
}

function plankDEGF(x) {
    return x *((-2 + x)/(-1 + x) + 2*Math.log(1.0 / (1.0 - x)));
}

function plankDDEGF(x) {
    return (2 - x*x + 2*(-1 + x)*(x-1)* Math.log(1 /(1 - x)))/((1 - x)*(1-x));
}

// Wa = Seq>=1(Pl)
function gammaWagon(x) {
    // sequence
    let p = plankEGF(x);
    let k = geometricTruncated(p, 1);
    if (k == -1) {
        return -1;
    }
    let seq = []
    let tmp;
    for (let i = 0; i < k; i++) {
        tmp = gammaPlanks(x);
        if (tmp == -1) {
            return -1;
        }
         seq.push(tmp);
    }
    // return {"wagon": seq};
    return seq;
}


// Pa = CYC(X) * CYC(Z)
function gammaPassenger(x) {
    let tmp = gammaCycEndpoint(x, endpointPerson, 1);
    if (tmp == -1) {
        return -1;
    }
    let tmp1 = gammaCycEndpoint(x, endpointPerson, 1);
    if (tmp1 == -1) {
        return -1;
    }
    return {"passenger": [tmp, tmp1]};
}

// SET(Pa)
function gammaSet(x, endpoint) {
    let p = cycEGF(x) * cycEGF(x);
    let k = poisson(p);
    let set = []
    let tmp;
    for (let i = 0; i < k; i++) {
        tmp = gammaPassenger(x);
        if (tmp == -1) {
            return -1;
        }
        set.push(tmp);
    }
    return set;
}

// x<0.682451
// -x^2*(log(-1/(x - 1)) + 1)/(x^2*(log(-1/(x - 1)) + 1) - 1)
function wagonEGF(x) {
    return plankEGF(x) * 1.0 / (1.0 - plankEGF(x));
}

function wagonDEGF(x) {
    return (2*x*(Math.log(-1/(x - 1)) + 1) - x*x/(x - 1))*x*x*(Math.log(-1/(x - 1)) + 1)/Math.pow((x*x*(Math.log(-1/(x - 1)) + 1) - 1),2) - 2*x*(Math.log(-1/(x - 1)) + 1)/(x*x*(Math.log(-1/(x - 1)) + 1) - 1) + x*x/((x*x*(Math.log(-1/(x - 1)) + 1) - 1)*(x - 1));
}

function wagonDDEGF(x) {
    return -2*Math.pow((2*x*(Math.log(-1/(x - 1)) + 1) - x*x/(x - 1)),2)*x*x*(Math.log(-1/(x - 1)) + 1)/Math.pow((x*x*(Math.log(-1/(x - 1)) + 1) - 1),3) - x*x*(4*x/(x - 1) - x*x/Math.pow(x-1, 2) - 2*Math.log(-1/(x - 1)) - 2)*(Math.log(-1/(x - 1)) + 1)/Math.pow((x*x*(Math.log(-1/(x - 1)) + 1) - 1), 2) + 4*(2*x*(Math.log(-1/(x - 1)) + 1) - x*x/(x - 1))*x*(Math.log(-1/(x - 1)) + 1)/Math.pow((x*x*(Math.log(-1/(x - 1)) + 1) - 1),2) - 2*(2*x*(Math.log(-1/(x - 1)) + 1) - x*x/(x - 1))*x*x/(Math.pow((x*x*(Math.log(-1/(x - 1)) + 1) - 1),2)*(x - 1)) - 2*(Math.log(-1/(x - 1)) + 1)/(x*x*(Math.log(-1/(x - 1)) + 1) - 1) + 4*x/((x*x*(Math.log(-1/(x - 1)) + 1) - 1)*(x - 1)) - x*x/((x*x*(Math.log(-1/(x - 1)) + 1) - 1)*Math.pow(x-1, 2));
}


//  e^(log(-1/(x - 1))^2)
function setPassengerEGF(x) {
    return Math.exp(cycEGF(x) * cycEGF(x));
}

function setPassengerDEGF(x) {
    return -2*Math.exp(Math.log(-1/(x - 1))*Math.log(-1/(x - 1)))*Math.log(-1/(x - 1))/(x - 1);
}

function setPassengerDDEGF(x) {
    return 4*Math.exp(Math.log(-1/(x - 1))*Math.log(-1/(x - 1)))*Math.log(-1/(x - 1)) *Math.log(-1/(x - 1)) /Math.pow((x - 1), 2) + 2*Math.exp(Math.log(-1/(x - 1)) * Math.log(-1/(x - 1)))*Math.log(-1/(x - 1))/Math.pow((x - 1),2) + 2*Math.exp(Math.log(-1/(x - 1))*Math.log(-1/(x - 1)))/((x - 1)*(x-1));
}

// 1/(x^2*(log(-1/(x - 1)) + 1)*e^(log(-1/(x - 1))^2)/(x^2*(log(-1/(x - 1)) + 1) - 1) + 1)
// Seq(Wa * Set(Pa)) = 1 / (1 - (Wa * Set(Pa)))
function gammaSeqWagonSetPassenger(x) {
    let p = wagonEGF(x) * setPassengerEGF(x);
    let k = geometric(p, 0);
    if (k == -1) {
        return -1;
    }
    let seq = []
    let tmp;
    let tmp1;
    for (let i = 0; i < k; i++) {
        tmp = gammaWagon(x);
        if (tmp == -1) {
            return -1;
        }
        tmp1 = gammaSet(x);
        if (tmp1 == -1) {
            return -1;
        }
        seq.push({"wagon": tmp, "passengers": tmp1});
    }
    return seq;
}

function trainEGF(x) {
    return -x*x*(Math.log(-1/(x - 1)) + 1)/((x*x*(Math.log(-1/(x - 1)) + 1) - 1)*(x*x*(Math.log(-1/(x - 1)) + 1)*Math.exp(Math.log(-1/(x - 1))*Math.log(-1/(x - 1)))/(x*x*(Math.log(-1/(x - 1)) + 1) - 1) + 1));
}

function trainDEGF(x) {
    return -(2*(Math.pow(Math.log(-1/(x - 1)), 3) + 2*Math.pow(Math.log(-1/(x - 1)), 2) + Math.log(-1/(x - 1)))*Math.pow(x, 4)*Math.exp(Math.pow(Math.log(-1/(x - 1)), 2)) - x*x*(2*Math.log(-1/(x - 1)) + 1) + 2*x*(Math.log(-1/(x - 1)) + 1))/((Math.pow(Math.log(-1/(x - 1)), 2) + 2*Math.log(-1/(x - 1)) + 1)*Math.pow(x, 5) - (Math.pow(Math.log(-1/(x - 1)), 2) + 2*Math.log(-1/(x - 1)) + 1)*Math.pow(x, 4) - 2*x*x*(Math.log(-1/(x - 1)) + 1) + 2*x*x*(Math.log(-1/(x - 1)) + 1) + ((Math.pow(Math.log(-1/(x - 1)), 2) + 2*Math.log(-1/(x - 1)) + 1)*Math.pow(x, 5) - (Math.pow(Math.log(-1/(x - 1)), 2) + 2*Math.log(-1/(x - 1)) + 1)*Math.pow(x, 4))*Math.exp(2*Math.pow(Math.log(-1/(x - 1)), 2)) + 2*((Math.pow(Math.log(-1/(x - 1)), 2) + 2*Math.log(-1/(x - 1)) + 1)*Math.pow(x, 5) - (Math.pow(Math.log(-1/(x - 1)), 2) + 2*Math.log(-1/(x - 1)) + 1)*Math.pow(x, 4) - x*x*(Math.log(-1/(x - 1)) + 1) + x*x*(Math.log(-1/(x - 1)) + 1))*Math.exp(Math.pow(Math.log(-1/(x - 1)), 2)) + x - 1);
}

function trainDDEGF(x) {
    return -(2*(2*Math.pow(Math.log(-1/(x - 1)), 5) + 5*Math.pow(Math.log(-1/(x - 1)), 4) + 2*Math.pow(Math.log(-1/(x - 1)), 3) - 4*Math.pow(Math.log(-1/(x - 1)), 2) - 4*Math.log(-1/(x - 1)) - 1)*Math.pow(x, 6)*Math.exp(2*Math.pow(Math.log(-1/(x - 1)), 2)) + (6*Math.pow(Math.log(-1/(x - 1)), 2) + 7*Math.log(-1/(x - 1)) + 3)*Math.pow(x, 4) - 4*(3*Math.pow(Math.log(-1/(x - 1)), 2) + 5*Math.log(-1/(x - 1)) + 2)*x*x*x + (6*Math.pow(Math.log(-1/(x - 1)), 2) + 14*Math.log(-1/(x - 1)) + 5)*x*x - (2*(2*Math.pow(Math.log(-1/(x - 1)), 5) + 7*Math.pow(Math.log(-1/(x - 1)), 4) + 10*Math.pow(Math.log(-1/(x - 1)), 3) + 8*Math.pow(Math.log(-1/(x - 1)), 2) + 4*Math.log(-1/(x - 1)) + 1)*Math.pow(x, 6) - (4*Math.pow(Math.log(-1/(x - 1)), 4) - 6*Math.pow(Math.log(-1/(x - 1)), 3) - 8*Math.pow(Math.log(-1/(x - 1)), 2) + 5*Math.log(-1/(x - 1)) + 5)*Math.pow(x, 4) - 4*(4*Math.pow(Math.log(-1/(x - 1)), 3) + 5*Math.pow(Math.log(-1/(x - 1)), 2) - Math.log(-1/(x - 1)) - 2)*x*x*x - 6*(Math.pow(Math.log(-1/(x - 1)), 2) + 2*Math.log(-1/(x - 1)) + 1)*x*x)*Math.exp(Math.pow(Math.log(-1/(x - 1)), 2)) - 4*x*Math.log(-1/(x - 1)) + 2*Math.log(-1/(x - 1)) + 2)/((Math.pow(Math.log(-1/(x - 1)), 3) + 3*Math.pow(Math.log(-1/(x - 1)), 2) + 3*Math.log(-1/(x - 1)) + 1)*Math.pow(x, 8) - 2*(Math.pow(Math.log(-1/(x - 1)), 3) + 3*Math.pow(Math.log(-1/(x - 1)), 2) + 3*Math.log(-1/(x - 1)) + 1)*Math.pow(x, 7) + (Math.pow(Math.log(-1/(x - 1)), 3) - 3*Math.log(-1/(x - 1)) - 2)*Math.pow(x, 6) + 6*(Math.pow(Math.log(-1/(x - 1)), 2) + 2*Math.log(-1/(x - 1)) + 1)*Math.pow(x, 5) - 3*(Math.pow(Math.log(-1/(x - 1)), 2) + Math.log(-1/(x - 1)))*Math.pow(x, 4) - 6*x*x*x*(Math.log(-1/(x - 1)) + 1) + x*x*(3*Math.log(-1/(x - 1)) + 2) + ((Math.pow(Math.log(-1/(x - 1)), 3) + 3*Math.pow(Math.log(-1/(x - 1)), 2) + 3*Math.log(-1/(x - 1)) + 1)*Math.pow(x, 8) - 2*(Math.pow(Math.log(-1/(x - 1)), 3) + 3*Math.pow(Math.log(-1/(x - 1)), 2) + 3*Math.log(-1/(x - 1)) + 1)*Math.pow(x, 7) + (Math.pow(Math.log(-1/(x - 1)), 3) + 3*Math.pow(Math.log(-1/(x - 1)), 2) + 3*Math.log(-1/(x - 1)) + 1)*Math.pow(x, 6))*Math.exp(3*Math.pow(Math.log(-1/(x - 1)), 2)) + 3*((Math.pow(Math.log(-1/(x - 1)), 3) + 3*Math.pow(Math.log(-1/(x - 1)), 2) + 3*Math.log(-1/(x - 1)) + 1)*Math.pow(x, 8) - 2*(Math.pow(Math.log(-1/(x - 1)), 3) + 3*Math.pow(Math.log(-1/(x - 1)), 2) + 3*Math.log(-1/(x - 1)) + 1)*Math.pow(x, 7) + (Math.pow(Math.log(-1/(x - 1)), 3) + 2*Math.pow(Math.log(-1/(x - 1)), 2) + Math.log(-1/(x - 1)))*Math.pow(x, 6) + 2*(Math.pow(Math.log(-1/(x - 1)), 2) + 2*Math.log(-1/(x - 1)) + 1)*Math.pow(x, 5) - (Math.pow(Math.log(-1/(x - 1)), 2) + 2*Math.log(-1/(x - 1)) + 1)*Math.pow(x, 4))*Math.exp(2*Math.pow(Math.log(-1/(x - 1)), 2)) + 3*((Math.pow(Math.log(-1/(x - 1)), 3) + 3*Math.pow(Math.log(-1/(x - 1)), 2) + 3*Math.log(-1/(x - 1)) + 1)*Math.pow(x, 8) - 2*(Math.pow(Math.log(-1/(x - 1)), 3) + 3*Math.pow(Math.log(-1/(x - 1)), 2) + 3*Math.log(-1/(x - 1)) + 1)*Math.pow(x, 7) + (Math.pow(Math.log(-1/(x - 1)), 3) + Math.pow(Math.log(-1/(x - 1)), 2) - Math.log(-1/(x - 1)) - 1)*Math.pow(x, 6) + 4*(Math.pow(Math.log(-1/(x - 1)), 2) + 2*Math.log(-1/(x - 1)) + 1)*Math.pow(x, 5) - (2*Math.pow(Math.log(-1/(x - 1)), 2) + 3*Math.log(-1/(x - 1)) + 1)*Math.pow(x, 4) - 2*x*x*x*(Math.log(-1/(x - 1)) + 1) + x*x*(Math.log(-1/(x - 1)) + 1))*Math.exp(Math.pow(Math.log(-1/(x - 1)), 2)) + 2*x - 1)
}
// dla tego zbieżność trzeba policzyć ekpreymentalnie - wolfram nie daje rady
// zbiega dla x <= 0.48512
// czyli kiedy wagonEGF(x) * setPassengerEGF(x) daje pbb <= 1
// Tr = Wa * SEQ(Wa * SET(Pa))
function gammaTrain(x) {
    let tmp = gammaWagon(x);
    if (tmp == -1) {
        return -1;
    }
    let tmp1 = gammaSeqWagonSetPassenger(x);
    if (tmp1 == -1 ) {
        return -1;
    }
    let train = [{"wagon": tmp, "passengers": []}];
    train = train.concat(tmp1);
    return {"wagons": train};
}


function test(fileName, rep, startX, endX, step) {
    let train;
    let stream = fs.createWriteStream("./results/" + fileName + ".txt", {flags:'a'});
    let streams = [];
    for (let i = 0; i < 8; i++) {
        streams.push(fs.createWriteStream("./results/" + fileName + "hist" + i + ".txt", {flags:'a'}));
    }
    let x = startX;
    let expected;
    let expVar;
    let startTime;
    let endTime;
    // statystyki: liczba wagonow, liczba plankow w wagonie, liczba kol,
    // endpointow w kole, liczba pasazerow w wagonie, liczba endpointow w glowie, liczba endpointow w brzuchu
    let stats = Array(8).fill(0.0);
    let min = Array(8).fill(0.0);
    let max = Array(8).fill(0.0);
    let variance = Array(8).fill(0.0);
    let res = Array(8).fill([]);
    let j = 0;
    let lastPrint = 0.0;
    let currLen;
    // tu jakas petla od x
    while (x <= endX) {
        for (let k = 0; k < stats.length; k++) {
            stats[k] = 0.0;
            min[k] = -1.0;
            max[k] = 0.0;
            variance[k] = 0.0;
            res[k] = [];
        }
        j = 0;
        if (x - lastPrint >= 0.01) {
            console.log(x);
            lastPrint = x;
        }
        while (j < rep) {
            startTime = performance?.timing ?? performance.timeOrigin;
            train = gammaTrain(x);
            endTime = (performance.now() + startTime) * 1000;
            if (train != -1) {
                currLen = 0;
                stats[7] += endTime;
                res[7].push(endTime);
                if (min[7] < 0 || min[7] > endTime) {
                    min[7] = endTime;
                }
                if (max[7] < endTime) {
                    max[7] = endTime;
                }

                /*stats[0] += train["wagons"].length;
                res[0].push(train["wagons"].length);
                if (min[0] < 0 || min[0] > train["wagons"].length) {
                    min[0] = train["wagons"].length;
                }
                if (max[0] < train["wagons"].length) {
                    max[0] = train["wagons"].length;
                }*/

                for (let k = 0; k < train["wagons"].length; k++) {
                    /*stats[1] += train["wagons"][k]["wagon"].length;
                    res[1].push(train["wagons"][k]["wagon"].length);
                    if (min[1] < 0 || min[1] > train["wagons"][k]["wagon"].length) {
                        min[1] = train["wagons"][k]["wagon"].length;
                    }
                    if (max[1] < train["wagons"][k]["wagon"].length) {
                        max[1] = train["wagons"][k]["wagon"].length;
                    }*/
                    for (let l = 0; l < train["wagons"][k]["wagon"].length; l++) {
                        /*stats[2] += train["wagons"][k]["wagon"][l]["wheels"].length;
                        res[2].push(train["wagons"][k]["wagon"][l]["wheels"].length);
                        if (min[2] < 0 || min[2] > train["wagons"][k]["wagon"][l]["wheels"].length) {
                            min[2] = train["wagons"][k]["wagon"][l]["wheels"].length;
                        }
                        if (max[2] < train["wagons"][k]["wagon"][l]["wheels"].length) {
                            max[2] = train["wagons"][k]["wagon"][l]["wheels"].length;
                        } */
                        for (let o = 0; o < train["wagons"][k]["wagon"][l]["wheels"].length; o++) {
                            currLen += train["wagons"][k]["wagon"][l]["wheels"][o].length;
                            /*stats[3] += train["wagons"][k]["wagon"][l]["wheels"][o].length;
                            res[3].push(train["wagons"][k]["wagon"][l]["wheels"][o].length);
                            if (min[3] < 0 || min[3] > train["wagons"][k]["wagon"][l]["wheels"][o].length) {
                                min[3] = train["wagons"][k]["wagon"][l]["wheels"][o].length;
                            }
                            if (max[3] < train["wagons"][k]["wagon"][l]["wheels"][o].length) {
                                max[3] = train["wagons"][k]["wagon"][l]["wheels"][o].length;
                            }*/
                        }
                    }

                    /*stats[4] += train["wagons"][k]["passengers"].length;
                    res[4].push(train["wagons"][k]["passengers"].length);
                    if (min[4] < 0 || min[4] > train["wagons"][k]["passengers"].length) {
                        min[4] = train["wagons"][k]["passengers"].length;
                    }
                    if (max[4] < train["wagons"][k]["passengers"].length) {
                        max[4] = train["wagons"][k]["passengers"].length;
                    } */
                    for (let l = 0; l < train["wagons"][k]["passengers"].length; l++) {
                        currLen += train["wagons"][k]["passengers"][l]["passenger"][0].length;
                        currLen += train["wagons"][k]["passengers"][l]["passenger"][1].length;
                        /*stats[5] += train["wagons"][k]["passengers"][l]["passenger"][0].length;
                        res[5].push(train["wagons"][k]["passengers"][l]["passenger"][0].length);
                        if (min[5] < 0 || min[5] > train["wagons"][k]["passengers"][l]["passenger"][0].length) {
                            min[5] = train["wagons"][k]["passengers"][l]["passenger"][0].length;
                        }
                        if (max[5] < train["wagons"][k]["passengers"][l]["passenger"][0].length) {
                            max[5] = train["wagons"][k]["passengers"][l]["passenger"][0].length;
                        }
                        stats[6] += train["wagons"][k]["passengers"][l]["passenger"][1].length;
                        res[6].push(train["wagons"][k]["passengers"][l]["passenger"][1].length);
                        if (min[6] < 0 || min[6] > train["wagons"][k]["passengers"][l]["passenger"][1].length) {
                            min[6] = train["wagons"][k]["passengers"][l]["passenger"][1].length;
                        }
                        if (max[6] < train["wagons"][k]["passengers"][l]["passenger"][1].length) {
                            max[6] = train["wagons"][k]["passengers"][l]["passenger"][1].length;
                        } */
                    }
                }
                stats[0] += currLen;
                res[0].push(currLen);
                if (min[0] < 0 || min[0] > currLen) {
                    min[0] = currLen;
                }
                if (max[0] < currLen) {
                    max[0] = currLen;
                }
                j += 1;
            } else {
                // powtórz
            }
        }
        // calc variance
        let avg = 0;
        for (let s = 0; s < stats.length; s++) {
            if (s == 0 || s == stats.length - 1) {
                avg = stats[s] / rep;
                for (let r = 0; r < res[s].length; r++) {
                    variance[s] += (res[s][r] - avg) * (res[s][r] - avg);
                    if (r != res[s].length - 1) {
                        streams[s].write(res[s][r] + " ");
                    } else {
                        streams[s].write(res[s][r].toString() + "\n");
                    }
                }
                variance[s] = variance[s] / (rep - 1);
            }
        }
        // calculate expected size from equation
        expected = x * trainDEGF(x) / trainEGF(x);
        expVar = ((x*x* trainDDEGF(x) / trainEGF(x)) + expected) - expected * expected;
        stream.write(x + " " + min[0] + " " + max[0] + " " + stats[0] / rep + " " + variance[0] + " " + expected + " " + expVar + " ");
        stream.write(min[stats.length - 1] + " " + max[stats.length - 1] + " " +
            stats[stats.length - 1] / rep + " " + variance[stats.length - 1] +  "\n");
        x += step;
    }
    stream.end();
    for (let i = 0; i < 8; i++) {
        streams[i].end();
    }
}
function testWagons(fileName, rep, startX, endX, step) {
    let wagon;
    let startTime;
    let endTime;
    let stream = fs.createWriteStream("./results/" + fileName + ".txt", {flags:'a'});
    let x = startX;
    let expected;
    let expVar;
    // statystyki: liczba plankow w wagonie, liczba kol, endpointow w kole
    let stats = Array(4).fill(0.0);
    let min = Array(4).fill(0.0);
    let max = Array(4).fill(0.0);
    let variance = Array(4).fill(0.0);
    let res = Array(4).fill([]);
    let streams = [];
    for (let i = 0; i < stats.length; i++) {
        streams.push(fs.createWriteStream("./results/" + fileName + "hist" + i + ".txt", {flags:'a'}));
    }
    let j = 0;
    let lastPrint = 0.0;
    let currLen;
    // tu jakas petla od x
    while (x <= endX) {
        for (let k = 0; k < stats.length; k++) {
            stats[k] = 0.0;
            min[k] = -1.0;
            max[k] = 0.0;
            variance[k] = 0.0;
            res[k] = [];
        }
        j = 0;
        if (x - lastPrint >= 0.1) {
            console.log(x);
            lastPrint = x;
        }
        while (j < rep) {
            startTime = performance?.timing ?? performance.timeOrigin;
            wagon = gammaWagon(x);
            endTime = (performance.now() + startTime) * 1000;
            if (wagon != -1) {
                currLen = 0;
                stats[3] += endTime;
                res[3].push(endTime);
                if (min[3] < 0 || min[3] > endTime) {
                    min[3] = endTime;
                }
                if (max[3] < endTime) {
                    max[3] = endTime;
                }
                /*stats[0] += wagon.length;
                res[0].push(wagon.length);
                if (min[0] < 0 || min[0] > wagon.length) {
                    min[0] = wagon.length;
                }
                if (max[0] < wagon.length) {
                    max[0] = wagon.length;
                } */
                for (let l = 0; l < wagon.length; l++) {
                    /*stats[1] += wagon[l]["wheels"].length;
                    res[1].push(wagon[l]["wheels"].length);
                    if (min[1] < 0 || min[1] > wagon[l]["wheels"].length) {
                        min[1] = wagon[l]["wheels"].length;
                    }
                    if (max[1] < wagon[l]["wheels"].length) {
                        max[1] = wagon[l]["wheels"].length;
                    }*/
                    for (let o = 0; o < wagon[l]["wheels"].length; o++) {
                        currLen += wagon[l]["wheels"][o].length;
                        /*stats[2] += wagon[l]["wheels"][o].length;
                        res[2].push(wagon[l]["wheels"][o].length);
                        if (min[2] < 0 || min[2] > wagon[l]["wheels"][o].length) {
                            min[2] = wagon[l]["wheels"][o].length;
                        }
                        if (max[2] < wagon[l]["wheels"][o].length) {
                            max[2] = wagon[l]["wheels"][o].length;
                        } */
                    }
                }
                stats[0] += currLen;
                res[0].push(currLen);
                if (min[0] < 0 || min[0] > currLen) {
                    min[0] = currLen;
                }
                if (max[0] < currLen) {
                    max[0] = currLen;
                }
                j += 1;
            } else {
                // powtórz
            }
        }
        // calc variance
        let avg = 0;
        for (let s = 0; s < stats.length; s++) {
            if (s == 0 || s == stats.length - 1) {
                avg = stats[s] / rep;
                for (let r = 0; r < res[s].length; r++) {
                    variance[s] += (res[s][r] - avg) * (res[s][r] - avg);
                    if (r != res[s].length - 1) {
                        streams[s].write(res[s][r] + " ");
                    } else {
                        streams[s].write(res[s][r].toString() + "\n");
                    }
                }
                variance[s] = variance[s] / (rep - 1);
            }
        }

        // wartość oczekiwana rozmiaru dla każdego
        expected = x * wagonDEGF(x) / wagonEGF(x);
        expVar = ((x*x*wagonDDEGF(x)) / wagonEGF(x) + expected) - expected*expected;
        stream.write(x + " " + min[0] + " " + max[0] + " " + stats[0] / rep + " " + variance[0] + " " + expected + " " + expVar + " ");
        // for (let s = 1; s < stats.length - 1; s++) {
        //     if (s == 1) {
        //         expected = x * plankDEGF(x) / plankEGF(x);
        //         expVar = ((x*x*plankDDEGF(x)) / plankEGF(x) + expected) - expected*expected;
        //     } else {
        //         // liczba endpointow w kole
        //         expected = x * cyc5DEGF(x) / cyc5EGF(x);
        //         expVar = ((x*x*cyc5DDEGF(x)) / cyc5EGF(x) + expected) - expected*expected;
        //     }
        //     stream.write(min[s] + " " + max[s] + " " + stats[s] / rep + " " + variance[s] + " " + expected + " " + expVar + " ");
        // }
        // czas
        stream.write(min[stats.length - 1] + " " + max[stats.length - 1] + " " +
            stats[stats.length - 1] / rep + " " + variance[stats.length - 1] + "\n");
        x += step;
    }
    stream.end();
    for (let i = 0; i < stats.length; i++) {
        streams[i].end();
    }
}

function testPlanks(fileName, rep, startX, endX, step) {
    let plank;
    let startTime;
    let endTime;
    let stream = fs.createWriteStream("./results/" + fileName + ".txt", {flags:'a'});
    let x = startX;
    let expVar;
    // statystyki: liczba kol (rozmiar planka),endpointow w kole
    let stats = Array(3).fill(0.0);
    let min = Array(3).fill(0.0);
    let max = Array(3).fill(0.0);
    let streams = [];
    for (let i = 0; i < stats.length; i++) {
        streams.push(fs.createWriteStream("./results/" + fileName + "hist" + i + ".txt", {flags:'a'}));
    }
    let variance = Array(3).fill(0.0);
    let res = Array(3).fill([]);
    let j = 0;
    let lastPrint = 0.0;
    let expected;
    let currLen;
    // tu jakas petla od x
    while (x <= endX) {
        for (let k = 0; k < stats.length; k++) {
            stats[k] = 0.0;
            min[k] = -1.0;
            max[k] = 0.0;
            variance[k] = 0.0;
            res[k] = [];
        }
        j = 0;
        if (x - lastPrint >= 0.1) {
            console.log(x);
            lastPrint = x;
        }
        while (j < rep) {
            startTime = performance?.timing ?? performance.timeOrigin;
            plank = gammaPlanks(x);
            endTime = (performance.now() + startTime) * 1000;
            if (plank != -1) {
                currLen = 0;
                stats[2] += endTime;
                res[2].push(endTime);
                if (min[2] < 0 || min[2] > endTime) {
                    min[2] = endTime;
                }
                if (max[2] < endTime) {
                    max[2] = endTime;
                }
                /*stats[0] += plank["wheels"].length;
                res[0].push(plank["wheels"].length);
                if (min[0] < 0 || min[0] > plank["wheels"].length) {
                    min[0] = plank["wheels"].length;
                }
                if (max[0] < plank["wheels"].length) {
                    max[0] = plank["wheels"].length;
                }*/
                for (let o = 0; o < plank["wheels"].length; o++) {
                    /*stats[1] += plank["wheels"][o].length;
                    res[1].push(plank["wheels"][o].length);
                     */
                    currLen += plank["wheels"][o].length;
                    /*if (min[1] < 0 || min[1] > plank["wheels"][o].length) {
                        min[1] = plank["wheels"][o].length;
                    }
                    if (max[1] < plank["wheels"][o].length) {
                        max[1] = plank["wheels"][o].length;
                    } */
                }
                stats[0] += currLen;
                res[0].push(currLen);
                if (min[0] < 0 || min[0] > currLen) {
                    min[0] = currLen;
                }
                if (max[0] < currLen) {
                    max[0] = currLen;
                }
                j += 1;
            } else {
                // powtórz
            }
        }
        // calc variance
        let avg = 0;
        for (let s = 0; s < stats.length; s++) {
            if (s == 0 || s == stats.length - 1) {
                avg = stats[s] / rep;
                for (let r = 0; r < res[s].length; r++) {
                    variance[s] += (res[s][r] - avg) * (res[s][r] - avg);
                    if (r != res[s].length - 1) {
                        streams[s].write(res[s][r] + " ");
                    } else {
                        streams[s].write(res[s][r].toString() + "\n");
                    }
                }
                variance[s] = variance[s] / (rep - 1);
            }
        }
        // calculate expected size from equation
        expected = x * plankDEGF(x) / plankEGF(x);
        expVar = (x*x*plankDDEGF(x)+ x*plankDEGF(x))/ plankEGF(x) - x*plankDEGF(x)/ plankEGF(x) * x*plankDEGF(x)/ plankEGF(x);
        stream.write(x + " " + min[0] + " " + max[0] + " " + stats[0] / rep + " " + variance[0] + " " + expected + " " + expVar + " ");
        // for (let s = 1; s < stats.length - 1; s++) {
        //     stream.write(min[s] + " " + max[s] + " " + stats[s] / rep + " " + variance[s] + " ");
        // }
        /*expected = x * cyc5DEGF(x) / cyc5EGF(x);
        expVar = (x * x * cyc5DDEGF(x) / cyc5EGF(x) + expected) - expected * expected;
        if (min[1] < 0) {
            min[1] = 0;
        }
        stream.write(min[1] + " " + max[1] + " " + stats[1] / rep + " " + variance[1] + " " + expected + " " + expVar + " "); */
        stream.write(min[stats.length - 1] + " " + max[stats.length - 1] + " " +
            stats[stats.length - 1] / rep + " " + variance[stats.length - 1] + "\n");
        x += step;
    }
    stream.end();
    for (let i = 0; i < stats.length; i++) {
        streams[i].end();
    }
}

function testCyc5(fileName, rep, startX, endX, step) {
    let cyc;
    let startTime;
    let endTime;
    let stream = fs.createWriteStream("./results/" + fileName + ".txt", {flags:'a'});
    let x = startX;
    let expected;
    let expVar;
    // statystyki endpointow w kole
    let stats = Array(2).fill(0.0);
    let min = Array(2).fill(0.0);
    let max = Array(2).fill(0.0);
    let variance = Array(2).fill(0.0);
    let res = Array(2).fill([]);
    let streams = [];
    for (let i = 0; i < stats.length; i++) {
        streams.push(fs.createWriteStream("./results/" + fileName + "hist" + i + ".txt", {flags:'a'}));
    }
    let j = 0;
    let lastPrint = 0.0;
    let currLen;
    // tu jakas petla od x
    while (x <= endX) {
        for (let k = 0; k < stats.length; k++) {
            stats[k] = 0.0;
            min[k] = -1.0;
            max[k] = 0.0;
            variance[k] = 0.0;
            res[k] = [];
        }
        j = 0;
        if (x - lastPrint >= 0.1) {
            console.log(x);
            lastPrint = x;
        }
        while (j < rep) {
            startTime = performance?.timing ?? performance.timeOrigin;
            cyc = gammaCycEndpoint(x, endpointPlank, 5);
            endTime = (performance.now() + startTime) * 1000;
            if (cyc != -1) {
                stats[1] += endTime;
                res[1].push(endTime);
                if (min[1] < 0 || min[1] > endTime) {
                    min[1] = endTime;
                }
                if (max[1] < endTime) {
                    max[1] = endTime;
                }
                stats[0] += cyc.length;
                res[0].push(cyc.length);
                if (min[0] < 0 || min[0] > cyc.length) {
                    min[0] = cyc.length;
                }
                if (max[0] < cyc.length) {
                    max[0] = cyc.length;
                }
                j += 1;
            } else {
                // powtórz
            }
        }
        // calc variance
        let avg = 0;
        for (let s = 0; s < stats.length; s++) {
            avg = stats[s] / rep;
            for (let r = 0; r < res[s].length; r++) {
                variance[s] += (res[s][r] - avg) * (res[s][r] - avg);
                if (r != res[s].length - 1) {
                    streams[s].write(res[s][r] + " ");
                } else {
                    streams[s].write(res[s][r].toString() + "\n");
                }
            }
            variance[s] = variance[s] / (rep - 1);
        }
        // calculate expected size from equation
        expected = x * cyc5DEGF(x) / cyc5EGF(x);
        expVar = (x * x * cyc5DDEGF(x) / cyc5EGF(x) + expected) - expected * expected;
        stream.write(x + " " + min[0] + " " + max[0] + " " +
            stats[0] / rep + " " + variance[0] + " " + expected + " " + expVar +" ");

        stream.write(min[stats.length - 1] + " " + max[stats.length - 1] + " " +
            stats[stats.length - 1] / rep + " " + variance[stats.length - 1] + " " + expected + " " + expVar +"\n");
        x += step;
    }
    stream.end();
    for (let i = 0; i < stats.length; i++) {
        streams[i].end();
    }
}

function testCyc2(fileName, rep, startX, endX, step) {
    let cyc;
    let startTime;
    let endTime;
    let stream = fs.createWriteStream("./results/" + fileName + ".txt", {flags:'a'});
    let x = startX;
    let expected;
    let expVar;
    // statystyki endpointow w kole
    let stats = Array(2).fill(0.0);
    let min = Array(2).fill(0.0);
    let max = Array(2).fill(0.0);
    let variance = Array(2).fill(0.0);
    let streams = [];
    for (let i = 0; i < stats.length; i++) {
        streams.push(fs.createWriteStream("./results/" + fileName + "hist" + i + ".txt", {flags:'a'}));
    }
    let res = Array(2).fill([]);
    let j = 0;
    let lastPrint = 0.0;
    // tu jakas petla od x
    while (x <= endX) {
        for (let k = 0; k < stats.length; k++) {
            stats[k] = 0.0;
            min[k] = -1.0;
            max[k] = 0.0;
            variance[k] = 0.0;
            res[k] = [];
        }
        j = 0;
        if (x - lastPrint >= 0.1) {
            console.log(x);
            lastPrint = x;
        }
        while (j < rep) {
            startTime = performance?.timing ?? performance.timeOrigin;
            cyc = gammaCycEndpoint(x, endpointPlank, 2);
            endTime = (performance.now() + startTime);
            if (cyc != -1) {
                stats[1] += endTime;
                res[1].push(endTime);
                if (min[1] < 0 || min[1] > endTime) {
                    min[1] = endTime;
                }
                if (max[1] < endTime) {
                    max[1] = endTime;
                }
                stats[0] += cyc.length;
                res[0].push(cyc.length);
                if (min[0] < 0 || min[0] > cyc.length) {
                    min[0] = cyc.length;
                }
                if (max[0] < cyc.length) {
                    max[0] = cyc.length;
                }
                j += 1;
            } else {
                // powtórz
            }
        }
        // calc variance
        let avg = 0;
        for (let s = 0; s < stats.length; s++) {
            avg = stats[s] / rep;
            for (let r = 0; r < res[s].length; r++) {
                variance[s] += (res[s][r] - avg) * (res[s][r] - avg);
                if (r != res[s].length - 1) {
                    streams[s].write(res[s][r] + " ");
                } else {
                    streams[s].write(res[s][r].toString() + "\n");
                }
            }
            variance[s] = variance[s] / (rep - 1);
        }
        // calculate expected size from equation
        expected = x * cyc2DEGF(x) / cyc2EGF(x);
        expVar = (x * x * cyc2DDEGF(x) / cyc2EGF(x) + expected) - expected * expected;
        stream.write(x + " " + min[0] + " " + max[0] + " " +
            stats[0] / rep + " " + variance[0] + " " + expected + " " + expVar+ " ");
        stream.write(min[stats.length - 1] + " " + max[stats.length - 1] + " " +
            stats[stats.length - 1] / rep + " " + variance[stats.length - 1] + "\n");
        x += step;
    }
    stream.end();
    for (let i = 0; i < stats.length; i++) {
        streams[i].end();
    }
}

function testCyc3(fileName, rep, startX, endX, step) {
    let cyc;
    let startTime;
    let endTime;
    let stream = fs.createWriteStream("./results/" + fileName + ".txt", {flags:'a'});
    let streams = [];
    for (let i = 0; i < 2; i++) {
        streams.push(fs.createWriteStream("./results/" + fileName + "hist" + i + ".txt", {flags:'a'}));
    }
    let x = startX;
    let expected;
    let expVar;
    // statystyki endpointow w kole
    let stats = Array(2).fill(0.0);
    let min = Array(2).fill(0.0);
    let max = Array(2).fill(0.0);
    let variance = Array(2).fill(0.0);
    let res = Array(2).fill([]);
    let j = 0;
    let lastPrint = 0.0;
    // tu jakas petla od x
    while (x <= endX) {
        for (let k = 0; k < stats.length; k++) {
            stats[k] = 0.0;
            min[k] = -1.0;
            max[k] = 0.0;
            variance[k] = 0.0;
            res[k] = [];
        }
        j = 0;
        if (x - lastPrint >= 0.1) {
            console.log(x);
            lastPrint = x;
        }
        while (j < rep) {
            startTime = performance?.timing ?? performance.timeOrigin;
            cyc = gammaCycEndpoint(x, endpointPlank, 3);
            endTime = (performance.now() + startTime) * 1000;
            if (cyc != -1) {
                stats[1] += endTime;
                res[1].push(endTime);
                if (min[1] < 0 || min[1] > endTime) {
                    min[1] = endTime;
                }
                if (max[1] < endTime) {
                    max[1] = endTime;
                }
                stats[0] += cyc.length;
                res[0].push(cyc.length);
                if (min[0] < 0 || min[0] > cyc.length) {
                    min[0] = cyc.length;
                }
                if (max[0] < cyc.length) {
                    max[0] = cyc.length;
                }
                j += 1;
            } else {
                // powtórz
            }
        }
        // calc variance
        let avg = 0;
        for (let s = 0; s < stats.length; s++) {
            avg = stats[s] / rep;
            for (let r = 0; r < res[s].length; r++) {
                variance[s] += (res[s][r] - avg) * (res[s][r] - avg);
                if (r != res[s].length - 1) {
                    streams[s].write(res[s][r] + " ");
                } else {
                    streams[s].write(res[s][r].toString() + "\n");
                }
            }
            variance[s] = variance[s] / (rep - 1);
        }
        // calculate expected size from equation
        expected = x * cyc3DEGF(x) / cyc3EGF(x);
        expVar = (x * x * cyc3DDEGF(x) / cyc3EGF(x) + expected) - expected * expected;
        stream.write(x + " " + min[0] + " " + max[0] + " " +
            stats[0] / rep + " " + variance[0] + " " + expected + " " + expVar + " ");
        stream.write(min[stats.length - 1] + " " + max[stats.length - 1] + " " +
            stats[stats.length - 1] / rep + " " + variance[stats.length - 1] + "\n");
        x += step;
    }
    stream.end();
    for (let i = 0; i < 2; i++) {
        streams[i].end();
    }
}

function testSetPassengers(fileName, rep, startX, endX, step) {
    let pass;
    let stream = fs.createWriteStream("./results/" + fileName + ".txt", {flags:'a'});
    let x = startX;
    let expected;
    let expVar;
    let startTime;
    let endTime;
    // statystyki: liczba pasazerow w wagonie, liczba endpointow w glowie, liczba endpointow w brzuchu
    let stats = Array(4).fill(0.0);
    let min = Array(4).fill(0.0);
    let max = Array(4).fill(0.0);
    let variance = Array(4).fill(0.0);
    let res = Array(4).fill([]);
    let streams = [];
    for (let i = 0; i < stats.length; i++) {
        streams.push(fs.createWriteStream("./results/" + fileName + "hist" + i + ".txt", {flags:'a'}));
    }
    let j = 0;
    let lastPrint = 0.0;
    let currLen;
    // tu jakas petla od x
    while (x <= endX) {
        for (let k = 0; k < stats.length; k++) {
            stats[k] = 0.0;
            min[k] = -1.0;
            max[k] = 0.0;
            variance[k] = 0.0;
            res[k] = [];
        }
        j = 0;
        if (x - lastPrint >= 0.1) {
            console.log(x);
            lastPrint = x;
        }
        while (j < rep) {
            startTime = performance?.timing ?? performance.timeOrigin;
            pass = gammaSet(x, endpointPerson);
            endTime = (performance.now() + startTime) * 1000;
            if (pass != -1) {
                currLen = 0;
                stats[3] += endTime;
                res[3].push(endTime);
                if (min[3] < 0 || min[3] > endTime) {
                    min[3] = endTime;
                }
                if (max[3] < endTime) {
                    max[3] = endTime;
                }
                /*stats[0] += pass.length;
                res[0].push(pass.length);
                if (min[0] < 0 || min[0] > pass.length) {
                    min[0] = pass.length;
                }
                if (max[0] < pass.length) {
                    max[0] = pass.length;
                } */
                for (let l = 0; l < pass.length; l++) {
                    currLen += pass[l]["passenger"][0].length;
                    currLen += pass[l]["passenger"][1].length;
                    /*stats[1] += pass[l]["passenger"][0].length;
                    res[1].push(pass[l]["passenger"][0].length);
                    if (min[1] < 0 || min[1] > pass[l]["passenger"][0].length) {
                        min[1] = pass[l]["passenger"][0].length;
                    }
                    if (max[1] < pass[l]["passenger"][0].length) {
                        max[1] = pass[l]["passenger"][0].length;
                    }
                    stats[2] += pass[l]["passenger"][1].length;
                    res[2].push(pass[l]["passenger"][1].length);
                    if (min[2] < 0 || min[2] > pass[l]["passenger"][1].length) {
                        min[2] = pass[l]["passenger"][1].length;
                    }
                    if (max[2] < pass[l]["passenger"][1].length) {
                        max[2] = pass[l]["passenger"][1].length;
                    } */
                }
                stats[0] += currLen;
                res[0].push(currLen);
                if (min[0] < 0 || min[0] > currLen) {
                    min[0] = currLen;
                }
                if (max[0] < currLen) {
                    max[0] = currLen;
                }
                j += 1;
            } else {
                // powtórz
            }
        }
        // calc variance
        let avg = 0;
        for (let s = 0; s < stats.length; s++) {
            if (s == 0 || s == stats.length - 1) {
                avg = stats[s] / rep;
                for (let r = 0; r < res[s].length; r++) {
                    variance[s] += (res[s][r] - avg) * (res[s][r] - avg);
                    if (r != res[s].length - 1) {
                        streams[s].write(res[s][r] + " ");
                    } else {
                        streams[s].write(res[s][r].toString() + "\n");
                    }
                }
                variance[s] = variance[s] / (rep - 1);
            }
        }
        // calculate expected size from equation
        expected = x * setPassengerDEGF(x) / setPassengerEGF(x);
        expVar = ((x*x* setPassengerDDEGF(x) / setPassengerEGF(x)) + expected) - expected * expected;
        stream.write(x + " " + min[0] + " " + max[0] + " " + stats[0] / rep + " " + variance[0] + " " + expected + " " + expVar + " ");
        // for (let s = 1; s < stats.length - 1; s++) {
        //     if (s == 1) {
        //         expected = x * cyc2DEGF(x) / cyc2EGF(x);
        //         expVar = ((x*x* cyc2DDEGF(x) / cyc2EGF(x)) + expected) - expected * expected;
        //     } else if (s == 2) {
        //         expected = x * cyc3DEGF(x) / cyc3EGF(x);
        //         expVar = ((x*x* cyc3DDEGF(x) / cyc3EGF(x)) + expected) - expected * expected;
        //     }
        //     stream.write(min[s] + " " + max[s] + " " + stats[s] / rep + " " + variance[s] + " " + expected + " " + expVar + " ");
        // }
        stream.write(min[stats.length - 1] + " " + max[stats.length - 1] + " " +
            stats[stats.length - 1] / rep + " " + variance[stats.length - 1] +"\n");
        x += step;
    }
    stream.end();
}



//console.log(JSON.stringify(gammaPlanks(0.8), null, " "))
// console.log(JSON.stringify(gammaWagon(0.01), null, " "));
let boundTrain =  0.48512;
let boundWagon = 0.682451;
// test("test_train", 100000, 0.57, 0.65, 0.001);
// testSetPassengers("passengers", 10000, 0.6, 0.9999, 0.001);
// console.log(gammaCycEndpoint(0.99999, endpointWheel, 3));
module.exports = {gammaTrain};
