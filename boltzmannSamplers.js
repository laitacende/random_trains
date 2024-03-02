const {loga, geometric, bernoulli, poisson, logaTruncated, geometricTruncated} = require("./distributions.js");
const fs = require('node:fs');
let endpointPlank = "○";
let endpointWheel = "•";
let endpointPerson = "p"

function endpointEGF(x) {
    return x; // rozmiar 1,  tworząca == 1 * x/1! - jeden element wielkości 1 i jedno możliwe wartościowanie
}

// function endpointGF(x) {
//     return x; // 1 * x
// }

function gammaEndpoint(x, symbol) {
    return symbol;
}


function gammaCycEndpoint(x, symbol, start) {
    let k = logaTruncated(endpointEGF(x), start);
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

function cyc5EGF(x) {
    return 1.0 / 12.0 * (-3.0 * Math.pow(x, 4) - 4.0 * x * x * x -6*x*x - 12*x -12*Math.log(1.0 - x));
}

function cyc5DEGF(x) {
    return Math.pow(x, 4)/(1 - x);
}

function cyc5DDEGF(x) {
    return ((4 - 3*x)*x*x*x)/((1 - x)*(1 - x))
}

// 1 + cyc >= 5
// start == 5
function gammaDisjointUnion(x, symbol, start) {
    let p = 1.0 / (1.0 + cyc5EGF(x)); // Uwaga: prawdopodobieństwo, że kół nie będzie jest bardzo duże...
                                // maleje wraz z wzrostem x, ale później taki x jest za duży dla zbieżności train...
    let tmp = bernoulli(p);
    if (tmp == 1) {
        return [];
    } else {
        let tmp1 = gammaCycEndpoint(x, symbol, start);
        if (tmp1 != -1) {
            return [tmp1];
        } else {
            return -1;
        }
    }
}

// seq(cyc) - koła NIEAKTUALNE
function gammaSeqCyc(x, symbol, start) {
    let p = Math.log(1 / (1 - endpointEGF(x)));
    let k = geometric(p, 0);
    if (k == -1) {
        return -1;
    }
    let seq = []
    let tmp;
    for (let i = 0; i < k; i++) {
        tmp = gammaCycEndpoint(x, symbol, 1);
        if (tmp == -1) {
            return -1;
        }
        seq.push(tmp);
    }
    return seq;
}

// Z * Z * (1 + Cyc>=5(Z))
function gammaPlanks(x) {
    // two endpoints
    let plank = [gammaEndpoint(x, endpointPlank),
            gammaEndpoint(x, endpointPlank)];
    let res = gammaDisjointUnion(x, endpointWheel, 5);
    if (res == -1) {
        return -1;
    }
    return {"planks": plank, "wheels": res};
}

// Z * Z * Seq(Cyc(Z)) NIEKATUALNE
function gammaPlanks2(x) {
    let w = gammaSeqCyc(x, endpointWheel, 1);
    if (w == -1) {
        return -1;
    }
    return {"planks": [gammaEndpoint(endpointEGF(x), endpointPlank),
            gammaEndpoint(endpointEGF(x), endpointPlank)],
            "wheels": w};
}

function plankEGF(x) {
    return x * x * (1 + cyc5EGF(x));
}

function plankDEGF(x) {
    return (x*(12 - 24*x + 6*x*x + 2*x*x*x + Math.pow(x, 4) + 9*Math.pow(x, 5)
        + 12*(-1 + x) *Math.log(1-x)))/(6*(1 - x));
}

function plankDDEGF(x) {
    return (12 - 36*x + 30*x*x - 4*x*x*x - Math.pow(x, 4) + 50*Math.pow(x, 5) -
        45*Math.pow(x, 6) - 12*(-1 + x)*(x - 1)*Math.log(1-x))/(6*(1 - x)*(1-x));
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

function cyc2EGF(x) {
    return -x - Math.log(1-x);
}

function cyc2DEGF(x) {
    return x/(1 - x);
}

function cyc2DDEGF(x) {
    return 1 / ((1 - x)*(1-x));
}

function cyc3EGF(x) {
    return 1.0 / 2.0 * (-x*x - 2 * x - 2 * Math.log(1-x));
}

function cyc3DEGF(x) {
    return x*x/(1 - x);
}

function cyc3DDEGF(x) {
    return -1 + 1 / ((1 - x)*(1-x))
}

// Pa = CYC>=2(X) * CYC>=3(Z)
function gammaPassenger(x) {
    let tmp = gammaCycEndpoint(x, endpointPerson, 2);
    if (tmp == -1) {
        return -1;
    }
    let tmp1 = gammaCycEndpoint(x, endpointPerson, 3);
    if (tmp1 == -1) {
        return -1;
    }
    return {"passenger": [tmp, tmp1]};
}

// SET(Pa)
function gammaSet(x, endpoint) {
    let p = cyc2EGF(x) * cyc3EGF(x);
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

// -1.09482<x<0.854909 zbieżne!
function wagonEGF(x) {
    return plankEGF(x) * 1.0 / (1.0 - plankEGF(x));
}

function wagonDEGF(x) {
    return (24*x* (12 - 24*x + 6*x*x + 2*x*x*x + Math.pow(x, 4) + 9*Math.pow(x, 5) + 12*(-1 + x)*Math.log(1 - x)))/((1 - x)*Math.pow((12 - 12*x*x + 12*x*x*x + 6*Math.pow(x, 4) + 4*Math.pow(x, 5) + 3*Math.pow(x, 6) + 12*x*x*Math.log(1 - x)),2));
}

function wagonDDEGF(x) {
    return (-48*x *(12 - 24*x + 6*x*x + 2*x*x*x + Math.pow(x, 4) + 9*Math.pow(x, 5) +
        12*(-1 + x)*Math.log(1 - x))*(-24*x + 36*x*x - (12*x*x)/(1 - x) + 24*x*x*x + 20*Math.pow(x, 4)
        + 18*Math.pow(x, 5) + 24*x*Math.log(1 - x)))/((1 - x)*Math.pow((12 - 12*x*x + 12*x*x*x + 6*Math.pow(x, 4)
        + 4*Math.pow(x, 5) + 3*Math.pow(x, 6) + 12*x*x*Math.log(1-x)), 3)) + (24*x* (-24 - (12* (-1 + x))/(1 - x) + 12*x
        + 6*x*x + 4*x*x*x + 45*Math.pow(x, 4) + 12*Math.log(1-x)))/((1 - x)* Math.pow((12 - 12*x*x + 12*x*x*x + 6*Math.pow(x, 4)
        + 4*Math.pow(x, 5) + 3*Math.pow(x, 6) + 12*x*x*Math.log(1-x)), 2)) + (24 *(12 - 24*x + 6*x*x +
        2*x*x*x + Math.pow(x, 4) + 9*Math.pow(x, 5) + 12*(-1 + x)* Math.log(1-x)))/((1 - x)* Math.pow((12 - 12*x*x
        + 12*x*x*x + 6*Math.pow(x, 4) + 4*Math.pow(x, 5) + 3*Math.pow(x, 6) + 12*x*x*Math.log(1-x)), 2)) +
    (24*x* (12 - 24*x + 6*x*x + 2*x*x*x + Math.pow(x, 4) + 9*Math.pow(x, 5) + 12*(-1 + x)*Math.log(1-x) ))/((1 - x)* (1-x) *Math.pow((12 -
        12*x*x + 12*x*x*x + 6*Math.pow(x, 4) + 4*Math.pow(x, 5) + 3*Math.pow(x, 6) + 12*x*x*Math.log(1-x)), 2));
}

function setPassengerEGF(x) {
    return Math.exp(cyc3EGF(x) * cyc2EGF(x));
}

function setPassengerDEGF(x) {
    return 1/2*(2*(x + 1/(x - 1) + 1)*(x + Math.log(-x + 1)) + (x*x + 2*x + 2*Math.log(-x + 1))*(1/(x - 1) + 1))*Math.exp(1/2*(x*x + 2*x + 2*Math.log(-x + 1))*(x + Math.log(-x + 1)));
}

function setPassengerDDEGF(x) {
    return 1/4*Math.pow((2*(x + 1/(x - 1) + 1)*(x + Math.log(-x + 1)) + (x*x + 2*x + 2*Math.log(-x + 1))*(1/(x - 1) + 1)), 2)*Math.exp(1/2*(x*x + 2*x + 2*Math.log(-x + 1))*(x + Math.log(-x + 1))) + 1/2*(4*(x + 1/(x - 1) + 1)*(1/(x - 1) + 1) - 2*(x + Math.log(-x + 1))*(1/((x-1)*(x-1)) - 1) - (x*x + 2*x + 2*Math.log(-x + 1))/((x-1)*(x-1)))*Math.exp(1/2*(x*x + 2*x + 2*Math.log(-x + 1))*(x + Math.log(-x + 1)))
}

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
    return (12*x*x - 12*x*x*x - 6*Math.pow(x, 4) - 4*Math.pow(x, 5) - 3*Math.pow(x, 6) - 12*x*x*Math.log(1 - x))/(12 - 12*x*x + 12*x*x*x + 6*Math.pow(x, 4) + 4*Math.pow(x, 5) + 
        3*Math.pow(x, 6) + 12*x*x*Math.log(1 - x)) * 1 / (1 - Math.exp(((-x - Math.log(1 - x)) *  (1/2*(-2*x - x*x - 2*Math.log(1 - x))))) * (12*x*x
        - 12*x*x*x - 6*Math.pow(x, 4) - 4*Math.pow(x, 5) - 3*Math.pow(x, 6) - 12 *x*x *Math.log(1 - x))/(12 - 12*x*x + 12*x*x*x + 6*Math.pow(x, 4) + 4*Math.pow(x, 5) + 3*Math.pow(x, 6) + 
        12 *x*x *Math.log(1 - x)));
}

function trainDEGF(x) {
    return 1/2*(3*Math.pow(x,6) + 4*Math.pow(x,5) + 6*Math.pow(x,4) + 12*x*x*x + 12*x*x*Math.log(-x + 1)
        - 12*x*x)*((3*Math.pow(x,6) + 4*Math.pow(x,5) + 6*Math.pow(x,4) + 12*x*x*x + 12*x*x*Math.log(-x + 1) - 12*x*x)*
        (2*(x + 1/(x - 1) + 1)*(x + Math.log(-x + 1)) + (x*x + 2*x + 2*Math.log(-x + 1))*(1/(x - 1) + 1))*Math.exp((1/2*(x*x + 2*x + 2*Math.log(-x + 1))
            *(x + Math.log(-x + 1))))/(3*Math.pow(x,6) + 4*Math.pow(x,5) + 6*Math.pow(x,4) + 12*x*x*x + 12*x*x*Math.log(-x + 1) - 12*x*x + 12)
        + 4*(9*Math.pow(x, 5) + 10*Math.pow(x, 4) + 12*x*x*x + 18*x*x + 12*x*Math.log(-x + 1) - 12*x + 6*x*x/(x - 1))*Math.exp((1/2*(x*x+ 2*x + 2*Math.log(-x + 1))*
            (x + Math.log(-x + 1))))/(3*Math.pow(x,6) + 4*Math.pow(x,5) + 6*Math.pow(x,4) + 12*x*x*x + 12*x*x*Math.log(-x + 1) -
            12*x*x + 12) - 4*(3*Math.pow(x,6) + 4*Math.pow(x,5) + 6*Math.pow(x,4) + 12*x*x*x + 12*x*x*Math.log(-x + 1) - 12*x*x)*
        (9*Math.pow(x,5) + 10*Math.pow(x,4) + 12*x*x*x + 18*x*x + 12*x*Math.log(-x + 1) - 12*x + 6*x*x/(x - 1))*
        Math.exp((1/2*(x*x + 2*x + 2*Math.log(-x + 1))*(x + Math.log(-x + 1))))/Math.pow((3*Math.pow(x,6) + 4*Math.pow(x,5) +
            6*Math.pow(x,4) + 12*x*x*x + 12*x*x*Math.log(-x + 1) - 12*x*x + 12), 2))/((3*Math.pow(x,6) + 4*Math.pow(x,5) +
        6*Math.pow(x,4) + 12*x*x*x + 12*x*x*Math.log(-x + 1) - 12*x*x + 12)*Math.pow(((3*Math.pow(x,6) + 4*Math.pow(x,5) + 6*Math.pow(x,4) +
        12*x*x*x + 12*x*x*Math.log(-x + 1) - 12*x*x)*Math.exp((1/2*(x*x + 2*x + 2*Math.log(-x + 1))*(x + Math.log(-x + 1))))/
        (3*Math.pow(x, 6) + 4*Math.pow(x, 5) + 6*Math.pow(x, 4) + 12*x*x*x + 12*x*x*Math.log(-x + 1) - 12*x*x + 12) + 1), 2))
        - 2*(9*Math.pow(x, 5) + 10*Math.pow(x, 4) + 12*x*x*x + 18*x*x + 12*x*Math.log(-x + 1) - 12*x + 6*x*x/(x - 1))/((3*Math.pow(x, 6)
            + 4*Math.pow(x, 5) + 6*Math.pow(x, 4) + 12*x*x*x + 12*x*x*Math.log(-x + 1) - 12*x*x + 12)*((3*Math.pow(x, 6)
            + 4*Math.pow(x, 5) + 6*Math.pow(x, 4) + 12*x*x*x + 12*x*x*Math.log(-x + 1) - 12*x*x)*
            Math.exp((1/2*(x*x + 2*x + 2*Math.log(-x + 1))*(x + Math.log(-x + 1))))/(3*Math.pow(x, 6) + 4*Math.pow(x, 5)
                + 6*Math.pow(x, 4) + 12*x*x*x + 12*x*x*Math.log(-x + 1) - 12*x*x + 12) + 1)) + 2*(3*Math.pow(x, 6) +
            4*Math.pow(x, 5) + 6*Math.pow(x, 4) + 12*x*x*x + 12*x*x*Math.log(-x + 1) - 12*x*x)*(9*Math.pow(x, 5) +
            10*Math.pow(x, 4) + 12*x*x*x + 18*x*x + 12*x*Math.log(-x + 1) - 12*x + 6*x*x/(x - 1))/(Math.pow((3*Math.pow(x, 6) +
            4*Math.pow(x, 5) + 6*Math.pow(x, 4) + 12*x*x*x + 12*x*x*Math.log(-x + 1) - 12*x*x + 12), 2)*((3*Math.pow(x, 6)
            + 4*Math.pow(x, 5) + 6*Math.pow(x, 4) + 12*x*x*x + 12*x*x*Math.log(-x + 1) - 12*x*x)*
            Math.exp((1/2*(x*x + 2*x + 2*Math.log(-x + 1))*(x + Math.log(-x + 1))))/(3*Math.pow(x, 6) + 4*Math.pow(x, 5)
                + 6*Math.pow(x, 4) + 12*x*x*x + 12*x*x*Math.log(-x + 1) - 12*x*x + 12) + 1));
}

function trainDDEGF(x) {
    return 1/4*(18144*Math.pow(x, 12) + 4032*Math.pow(x, 11) + 7200*Math.pow(x, 10) + 18432*Math.pow(x, 9) + 2304*Math.pow(x, 8)*(12*Math.log(-x + 1) - 19) - 13824*Math.pow(x, 7)*(Math.log(-x + 1) - 1) + 6912*Math.pow(x, 6)*(Math.log(-x + 1) - 9) + 2304*Math.pow(x, 5)*(12*Math.log(-x + 1) - 5) + 1152*(36*Math.pow(Math.log(-x + 1), 2) - 180*Math.log(-x + 1) + 179)*Math.pow(x, 4) - 4608*(18*Math.pow(Math.log(-x + 1), 2) - 54*Math.log(-x + 1) + 37)*x*x*x + 6912*(6*Math.pow(Math.log(-x + 1), 2) - 14*Math.log(-x + 1) + 11)*x*x + 13824*x*(2*Math.log(-x + 1) - 3) - (243*Math.pow(x, 24) + 324*Math.pow(x, 23)*(Math.log(-x + 1) + 4) + 54*(2*Math.pow(Math.log(-x + 1), 2) + 34*Math.log(-x + 1) + 77)*Math.pow(x, 22) + 72*(9*Math.pow(Math.log(-x + 1), 2) + 84*Math.log(-x + 1) + 155)*Math.pow(x, 21) + 6*(366*Math.pow(Math.log(-x + 1), 2) + 3272*Math.log(-x + 1) + 3401)*Math.pow(x, 20) + 4*(2512*Math.pow(Math.log(-x + 1), 2) + 11030*Math.log(-x + 1) + 6895)*Math.pow(x, 19) + 4*(324*Math.pow(Math.log(-x + 1), 3) + 7310*Math.pow(Math.log(-x + 1), 2) + 18785*Math.log(-x + 1) + 6963)*Math.pow(x, 18) + 16*(378*Math.pow(Math.log(-x + 1), 3) + 3787*Math.pow(Math.log(-x + 1), 2) + 7184*Math.log(-x + 1) + 382)*Math.pow(x, 17) + 16*(981*Math.pow(Math.log(-x + 1), 3) + 8010*Math.pow(Math.log(-x + 1), 2) + 6134*Math.log(-x + 1) - 1271)*Math.pow(x, 16) + 32*(1602*Math.pow(Math.log(-x + 1), 3) + 5283*Math.pow(Math.log(-x + 1), 2) + 863*Math.log(-x + 1) - 891)*Math.pow(x, 15) + 144*(36*Math.pow(Math.log(-x + 1), 4) + 676*Math.pow(Math.log(-x + 1), 3) + 973*Math.pow(Math.log(-x + 1), 2) - 264*Math.log(-x + 1) - 87)*Math.pow(x, 14) + 576*(30*Math.pow(Math.log(-x + 1), 4) + 204*Math.pow(Math.log(-x + 1), 3) + 179*Math.pow(Math.log(-x + 1), 2) - 195*Math.log(-x + 1) + 20)*Math.pow(x, 13) + 288*(102*Math.pow(Math.log(-x + 1), 4) + 558*Math.pow(Math.log(-x + 1), 3) - 301*Math.pow(Math.log(-x + 1), 2) - 224*Math.log(-x + 1) + 150)*Math.pow(x, 12) + 1728*(40*Math.pow(Math.log(-x + 1), 4) + 30*Math.pow(Math.log(-x + 1), 3) - 101*Math.pow(Math.log(-x + 1), 2) + 60*Math.log(-x + 1) + 6)*Math.pow(x, 11) + 1728*(4*Math.pow(Math.log(-x + 1), 5) + 38*Math.pow(Math.log(-x + 1), 4) - 63*Math.pow(Math.log(-x + 1), 3) + 10*Math.pow(Math.log(-x + 1), 2) + 88*Math.log(-x + 1) - 59)*Math.pow(x, 10) + 6912*(2*Math.pow(Math.log(-x + 1), 5) - Math.pow(Math.log(-x + 1), 4) - 8*Math.pow(Math.log(-x + 1), 3) + 32*Math.pow(Math.log(-x + 1), 2) - 33*Math.log(-x + 1) + 9)*Math.pow(x, 9) + 3456*(2*Math.pow(Math.log(-x + 1), 5) - 8*Math.pow(Math.log(-x + 1), 4) + 32*Math.pow(Math.log(-x + 1), 3) - 41*Math.pow(Math.log(-x + 1), 2) + 14*Math.log(-x + 1) + 1)*Math.pow(x, 8) + 6912*(2*Math.pow(Math.log(-x + 1), 4) - 2*Math.pow(Math.log(-x + 1), 3) - 3*Math.pow(Math.log(-x + 1), 2) + 4*Math.log(-x + 1) - 1)*Math.pow(x, 7) + 6912*(Math.pow(Math.log(-x + 1), 4) - 3*Math.pow(Math.log(-x + 1), 3) + 3*Math.pow(Math.log(-x + 1), 2) - Math.log(-x + 1))*Math.pow(x, 6))*Math.exp(x*x*x + x*x*Math.log(-x + 1) + 2*x*x + 4*x*Math.log(-x + 1) + 2*Math.pow(Math.log(-x + 1), 2)) + (243*Math.pow(x, 24) + 324*Math.pow(x, 23)*(Math.log(-x + 1) + 4) + 54*(2*Math.pow(Math.log(-x + 1), 2) + 34*Math.log(-x + 1) + 77)*Math.pow(x, 22) + 72*(9*Math.pow(Math.log(-x + 1), 2) + 84*Math.log(-x + 1) + 164)*Math.pow(x, 21) + 6*(366*Math.pow(Math.log(-x + 1), 2) + 3308*Math.log(-x + 1) + 3743)*Math.pow(x, 20) + 4*(2512*Math.pow(Math.log(-x + 1), 2) + 11138*Math.log(-x + 1) + 8137)*Math.pow(x, 19) + 4*(324*Math.pow(Math.log(-x + 1), 3) + 7310*Math.pow(Math.log(-x + 1), 2) + 18911*Math.log(-x + 1) + 10380)*Math.pow(x, 18) + 16*(378*Math.pow(Math.log(-x + 1), 3) + 3787*Math.pow(Math.log(-x + 1), 2) + 7801*Math.log(-x + 1) + 1175)*Math.pow(x, 17) + 16*(981*Math.pow(Math.log(-x + 1), 3) + 8199*Math.pow(Math.log(-x + 1), 2) + 6979*Math.log(-x + 1) - 367)*Math.pow(x, 16) + 16*(3204*Math.pow(Math.log(-x + 1), 3) + 10800*Math.pow(Math.log(-x + 1), 2) + 3332*Math.log(-x + 1) + 297)*Math.pow(x, 15) + 24*(216*Math.pow(Math.log(-x + 1), 4) + 4056*Math.pow(Math.log(-x + 1), 3) + 5996*Math.pow(Math.log(-x + 1), 2) + 2110*Math.log(-x + 1) - 877)*Math.pow(x, 14) + 48*(360*Math.pow(Math.log(-x + 1), 4) + 2448*Math.pow(Math.log(-x + 1), 3) + 3212*Math.pow(Math.log(-x + 1), 2) - 2434*Math.log(-x + 1) + 1695)*Math.pow(x, 13) + 48*(612*Math.pow(Math.log(-x + 1), 4) + 3636*Math.pow(Math.log(-x + 1), 3) - 1874*Math.pow(Math.log(-x + 1), 2) + 147*Math.log(-x + 1) + 3592)*Math.pow(x, 12) + 192*(360*Math.pow(Math.log(-x + 1), 4) + 294*Math.pow(Math.log(-x + 1), 3) - 771*Math.pow(Math.log(-x + 1), 2) + 1600*Math.log(-x + 1) - 342)*Math.pow(x, 11) + 576*(12*Math.pow(Math.log(-x + 1), 5) + 114*Math.pow(Math.log(-x + 1), 4) - 185*Math.pow(Math.log(-x + 1), 3) + 344*Math.pow(Math.log(-x + 1), 2) - 87*Math.log(-x + 1) + 10)*Math.pow(x, 10) + 576*(24*Math.pow(Math.log(-x + 1), 5) - 12*Math.pow(Math.log(-x + 1), 4) + 44*Math.pow(Math.log(-x + 1), 3) - 20*Math.pow(Math.log(-x + 1), 2) + 178*Math.log(-x + 1) - 135)*Math.pow(x, 9) + 576*(12*Math.pow(Math.log(-x + 1), 5) - 12*Math.pow(Math.log(-x + 1), 4) - 24*Math.pow(Math.log(-x + 1), 3) + 164*Math.pow(Math.log(-x + 1), 2) + 109*Math.log(-x + 1) - 466)*Math.pow(x, 8) + 6912*(2*Math.pow(Math.log(-x + 1), 3) + 24*Math.pow(Math.log(-x + 1), 2) - 77*Math.log(-x + 1) + 31)*Math.pow(x, 7) + 3456*(20*Math.pow(Math.log(-x + 1), 3) - 79*Math.pow(Math.log(-x + 1), 2) + 20*Math.log(-x + 1) + 20)*Math.pow(x, 6) - 6912*(2*Math.pow(Math.log(-x + 1), 3) + 23*Math.pow(Math.log(-x + 1), 2) - 38*Math.log(-x + 1) + 19)*Math.pow(x, 5) - 20736*(3*Math.pow(Math.log(-x + 1), 3) - 8*Math.pow(Math.log(-x + 1), 2) + 13*Math.log(-x + 1) - 10)*Math.pow(x, 4) - 82944*(Math.pow(Math.log(-x + 1), 2) - 3*Math.log(-x + 1) + 2)*x*x*x + 41472*(Math.pow(Math.log(-x + 1), 2) - 2*Math.log(-x + 1) + 1)*x*x)*Math.exp(1/2*x*x*x + 1/2*x*x*Math.log(-x + 1) + x*x + 2*x*Math.log(-x + 1) + Math.pow(Math.log(-x + 1), 2)) - 13824*Math.log(-x + 1) + 13824)/(27*Math.pow(x, 20) + 54*Math.pow(x, 19) + 117*Math.pow(x, 18) + 316*Math.pow(x, 17) + 2*Math.pow(x, 16)*(162*Math.log(-x + 1) - 91) + 4*Math.pow(x, 15)*(54*Math.log(-x + 1) - 11) + 36*Math.pow(x, 14)*(13*Math.log(-x + 1) - 1) + 72*Math.pow(x, 13)*(20*Math.log(-x + 1) - 33) + 36*(36*Math.pow(Math.log(-x + 1), 2) - 128*Math.log(-x + 1) + 99)*Math.pow(x, 12) - 144*(6*Math.pow(Math.log(-x + 1), 2) - 12*Math.log(-x + 1) - 1)*Math.pow(x, 11) + 144*(3*Math.pow(Math.log(-x + 1), 2) + 9*Math.log(-x + 1) - 32)*Math.pow(x, 10) + 1728*(Math.pow(Math.log(-x + 1), 2) - 6*Math.log(-x + 1) + 6)*Math.pow(x, 9) + 864*(2*Math.pow(Math.log(-x + 1), 3) - 15*Math.pow(Math.log(-x + 1), 2) + 31*Math.log(-x + 1) - 17)*Math.pow(x, 8) - 864*(4*Math.pow(Math.log(-x + 1), 3) - 18*Math.pow(Math.log(-x + 1), 2) + 20*Math.log(-x + 1) + 1)*Math.pow(x, 7) + 432*(4*Math.pow(Math.log(-x + 1), 3) - 48*Math.log(-x + 1) + 57)*Math.pow(x, 6) - 1728*(6*Math.pow(Math.log(-x + 1), 2) - 18*Math.log(-x + 1) + 11)*Math.pow(x, 5) + 2592*(2*Math.pow(Math.log(-x + 1), 2) - 2*Math.log(-x + 1) - 3)*Math.pow(x, 4) - 5184*x*x*x*(2*Math.log(-x + 1) - 3) + 1728*x*x*(3*Math.log(-x + 1) - 2) + (27*Math.pow(x, 20) + 54*Math.pow(x, 19) + 117*Math.pow(x, 18) + 316*Math.pow(x, 17) + 2*Math.pow(x, 16)*(162*Math.log(-x + 1) - 91) + 4*Math.pow(x, 15)*(54*Math.log(-x + 1) - 11) + 36*Math.pow(x, 14)*(13*Math.log(-x + 1) - 10) + 288*Math.pow(x, 13)*(5*Math.log(-x + 1) - 9) + 72*(18*Math.pow(Math.log(-x + 1), 2) - 64*Math.log(-x + 1) + 43)*Math.pow(x, 12) - 432*(2*Math.pow(Math.log(-x + 1), 2) - 4*Math.log(-x + 1) + 3)*Math.pow(x, 11) + 432*(Math.pow(Math.log(-x + 1), 2) - 3*Math.log(-x + 1))*Math.pow(x, 10) + 1728*(Math.pow(Math.log(-x + 1), 2) - 5*Math.log(-x + 1) + 5)*Math.pow(x, 9) + 864*(2*Math.pow(Math.log(-x + 1), 3) - 15*Math.pow(Math.log(-x + 1), 2) + 30*Math.log(-x + 1) - 17)*Math.pow(x, 8) - 1728*(2*Math.pow(Math.log(-x + 1), 3) - 9*Math.pow(Math.log(-x + 1), 2) + 12*Math.log(-x + 1) - 5)*Math.pow(x, 7) + 1728*(Math.pow(Math.log(-x + 1), 3) - 3*Math.pow(Math.log(-x + 1), 2) + 3*Math.log(-x + 1) - 1)*Math.pow(x, 6))*Math.exp(3/2*x*x*x + 3/2*x*x*Math.log(-x + 1) + 3*x*x + 6*x*Math.log(-x + 1) + 3*Math.pow(Math.log(-x + 1), 2)) + 3*(27*Math.pow(x, 20) + 54*Math.pow(x, 19) + 117*Math.pow(x, 18) + 316*Math.pow(x, 17) + 2*Math.pow(x, 16)*(162*Math.log(-x + 1) - 91) + 4*Math.pow(x, 15)*(54*Math.log(-x + 1) - 11) + 36*Math.pow(x, 14)*(13*Math.log(-x + 1) - 7) + 360*Math.pow(x, 13)*(4*Math.log(-x + 1) - 7) + 12*(108*Math.pow(Math.log(-x + 1), 2) - 384*Math.log(-x + 1) + 271)*Math.pow(x, 12) - 48*(18*Math.pow(Math.log(-x + 1), 2) - 36*Math.log(-x + 1) + 17)*Math.pow(x, 11) + 48*(9*Math.pow(Math.log(-x + 1), 2) - 9*Math.log(-x + 1) - 32)*Math.pow(x, 10) + 576*(3*Math.pow(Math.log(-x + 1), 2) - 16*Math.log(-x + 1) + 16)*Math.pow(x, 9) + 144*(12*Math.pow(Math.log(-x + 1), 3) - 90*Math.pow(Math.log(-x + 1), 2) + 182*Math.log(-x + 1) - 105)*Math.pow(x, 8) - 576*(6*Math.pow(Math.log(-x + 1), 3) - 27*Math.pow(Math.log(-x + 1), 2) + 34*Math.log(-x + 1) - 10)*Math.pow(x, 7) + 1728*(Math.pow(Math.log(-x + 1), 3) - 2*Math.pow(Math.log(-x + 1), 2) - 2*Math.log(-x + 1) + 4)*Math.pow(x, 6) - 3456*(Math.pow(Math.log(-x + 1), 2) - 3*Math.log(-x + 1) + 2)*Math.pow(x, 5) + 1728*(Math.pow(Math.log(-x + 1), 2) - 2*Math.log(-x + 1) + 1)*Math.pow(x, 4))*Math.exp(x*x*x + x*x*Math.log(-x + 1) + 2*x*x + 4*x*Math.log(-x + 1) + 2*Math.pow(Math.log(-x + 1), 2)) + 3*(27*Math.pow(x, 20) + 54*Math.pow(x, 19) + 117*Math.pow(x, 18) + 316*Math.pow(x, 17) + 2*Math.pow(x, 16)*(162*Math.log(-x + 1) - 91) + 4*Math.pow(x, 15)*(54*Math.log(-x + 1) - 11) + 36*Math.pow(x, 14)*(13*Math.log(-x + 1) - 4) + 144*Math.pow(x, 13)*(10*Math.log(-x + 1) - 17) + 48*(27*Math.pow(Math.log(-x + 1), 2) - 96*Math.log(-x + 1) + 71)*Math.pow(x, 12) - 48*(18*Math.pow(Math.log(-x + 1), 2) - 36*Math.log(-x + 1) + 7)*Math.pow(x, 11) + 48*(9*Math.pow(Math.log(-x + 1), 2) + 9*Math.log(-x + 1) - 64)*Math.pow(x, 10) + 576*(3*Math.pow(Math.log(-x + 1), 2) - 17*Math.log(-x + 1) + 17)*Math.pow(x, 9) + 144*(12*Math.pow(Math.log(-x + 1), 3) - 90*Math.pow(Math.log(-x + 1), 2) + 184*Math.log(-x + 1) - 105)*Math.pow(x, 8) - 288*(12*Math.pow(Math.log(-x + 1), 3) - 54*Math.pow(Math.log(-x + 1), 2) + 64*Math.log(-x + 1) - 9)*Math.pow(x, 7) + 144*(12*Math.pow(Math.log(-x + 1), 3) - 12*Math.pow(Math.log(-x + 1), 2) - 84*Math.log(-x + 1) + 109)*Math.pow(x, 6) - 576*(12*Math.pow(Math.log(-x + 1), 2) - 36*Math.log(-x + 1) + 23)*Math.pow(x, 5) + 864*(4*Math.pow(Math.log(-x + 1), 2) - 6*Math.log(-x + 1) - 1)*Math.pow(x, 4) - 1728*x*x*x*(2*Math.log(-x + 1) - 3) + 1728*x*x*(Math.log(-x + 1) - 1))*Math.exp(1/2*x*x*x + 1/2*x*x*Math.log(-x + 1) + x*x + 2*x*Math.log(-x + 1) + Math.pow(Math.log(-x + 1), 2)) - 3456*x + 1728);
}
// dla tego zbieżność trzeba policzyć ekpreymentalnie - wolfram nie daje rady
// zbiega dla x <= 0.6695314
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
let boundTrain = 0.6695314;
let boundWagon = 0.854909;
// test("test_train", 100000, 0.57, 0.65, 0.001);
// test("train", 10000, 0.3, boundTrain,0.001);
// gammaTrain(boundTrain);
module.exports = {gammaTrain};
// Z(x) = x
// Z'(x) = x
// CYC(Z)(x) = log(1 / (1 - Z'(x)) = log(1 / (1 - x)))
// SEQ(Z)(x) = 1 / (1 - Z'(x)) = 1 / (1 - x)
// SEQ(CYC(Z))(x) = 1 / (1 - log(1 / (1 - x)))

// Cyc>=5(x) = 1/12 (-12 x - 6 x^2 - 4 x^3 - 3 x^4 - 12 log(1 - x))
// Cyc>=5'(x) = x^4/(1 - x)
// Cyc>=5''(x) = ((4 - 3 x) x^3)/(1 - x)^2

// Cyc>=2(x) = -x - log(1 - x)
// Cyc>=2'(x) = x/(1 - x)
// Cyc>=2''(x) = (1 - x)^(-2)

// Cyc>=3(x) = 1/2 (-2 x - x^2 - 2 log(1 - x))
// Cyc>=3'(x) = x^2/(1 - x)
// Cyc>3''(x) = -1 + (1 - x)^(-2)

// Pl(x) = (1 + 1/12 (-12 x - 6 x^2 - 4 x^3 - 3 x^4 - 12 log(1 - x)))*x^2)
// Pl'(x) = (x (12 - 24 x + 6 x^2 + 2 x^3 + x^4 + 9 x^5 + 12 (-1 + x) Log[1 - x]))/(6 (1 - x))
// Pl''(x) = (12 - 36 x + 30 x^2 - 4 x^3 - x^4 + 50 x^5 - 45 x^6 - 12 (-1 + x)^2 Log[1 - x])/(6 (1 - x)^2)

// Wa(x) = (12*x^2 - 12*x^3 - 6*x^4 - 4*x^5 - 3*x^6 - 12*x^2*log(1 - x))/(12 - 12*x^2 + 12*x^3 + 6*x^4 + 4*x^5 + 3*x^6 + 12*x^2*log(1 - x))
// Wa'(x) = (24 x (12 - 24 x + 6 x^2 + 2 x^3 + x^4 + 9 x^5 + 12 (-1 + x) Log[1 - x]))/((1 - x) (12 - 12 x^2 + 12 x^3 + 6 x^4 + 4 x^5 + 3 x^6 + 12 x^2 Log[1 - x])^2)
// Wa''(x) = (-48 x (12 - 24 x + 6 x^2 + 2 x^3 + x^4 + 9 x^5 + 12 (-1 + x) Log[1 - x]) (-24 x + 36 x^2 - (12 x^2)/(1 - x) + 24 x^3 + 20 x^4 + 18 x^5 + 24 x Log[1 - x]))/((1 - x) (12 - 12 x^2 + 12 x^3 + 6 x^4 + 4 x^5 + 3 x^6 + 12 x^2 Log[1 - x])^3) + (24 x (-24 - (12 (-1 + x))/(1 - x) + 12 x + 6 x^2 + 4 x^3 + 45 x^4 + 12 Log[1 - x]))/((1 - x) (12 - 12 x^2 + 12 x^3 + 6 x^4 + 4 x^5 + 3 x^6 + 12 x^2 Log[1 - x])^2) + (24 (12 - 24 x + 6 x^2 + 2 x^3 + x^4 + 9 x^5 + 12 (-1 + x) Log[1 - x]))/((1 - x) (12 - 12 x^2 + 12 x^3 + 6 x^4 + 4 x^5 + 3 x^6 + 12 x^2 Log[1 - x])^2) + (24 x (12 - 24 x + 6 x^2 + 2 x^3 + x^4 + 9 x^5 + 12 (-1 + x) Log[1 - x]))/((1 - x)^2 (12 - 12 x^2 + 12 x^3 + 6 x^4 + 4 x^5 + 3 x^6 + 12 x^2 Log[1 - x])^2)

// Pa(x) = (-x - log(1 - x)) *  (1/2 (-2 x - x^2 - 2 log(1 - x)))

// SET(Pa)(x) = e^((-x - log(1 - x)) *  (1/2* (-2*x - x^2 - 2*log(1 - x))))
//  SET(Pa)'(x) = 1/2*(2*(x + 1/(x - 1) + 1)*(x + log(-x + 1)) + (x^2 + 2*x + 2*log(-x + 1))*(1/(x - 1) + 1))*e^(1/2*(x^2 + 2*x + 2*log(-x + 1))*(x + log(-x + 1)))
// SET(Pa)''(x) = 1/4*(2*(x + 1/(x - 1) + 1)*(x + log(-x + 1)) + (x^2 + 2*x + 2*log(-x + 1))*(1/(x - 1) + 1))^2*e^(1/2*(x^2 + 2*x + 2*log(-x + 1))*(x + log(-x + 1))) + 1/2*(4*(x + 1/(x - 1) + 1)*(1/(x - 1) + 1) - 2*(x + log(-x + 1))*(1/(x - 1)^2 - 1) - (x^2 + 2*x + 2*log(-x + 1))/(x - 1)^2)*e^(1/2*(x^2 + 2*x + 2*log(-x + 1))*(x + log(-x + 1)))

// Wa * set(Pa) = e^((-x - log(1 - x)) *  (1/2 (-2 x - x^2 - 2 log(1 - x)))) * (12 x^2 - 12 x^3 - 6 x^4 - 4 x^5 - 3 x^6 - 12 x^2 log(1 - x))/(12 - 12 x^2 + 12 x^3 + 6 x^4 + 4 x^5 + 3 x^6 + 12 x^2 log(1 - x))
// SEQ(WA*SET(PA)) = 1 / (1 - e^((-x - log(1 - x)) *  (1/2 (-2 x - x^2 - 2 log(1 - x)))) * (12 x^2 - 12 x^3 - 6 x^4 - 4 x^5 - 3 x^6 - 12 x^2 log(1 - x))/(12 - 12 x^2 + 12 x^3 + 6 x^4 + 4 x^5 + 3 x^6 + 12 x^2 log(1 - x)))

// Tr(x) = (12 x^2 - 12 x^3 - 6 x^4 - 4 x^5 - 3 x^6 - 12 x^2 log(1 - x))/(12 - 12 x^2 + 12 x^3 + 6 x^4 + 4 x^5 + 3 x^6 + 12 x^2 log(1 - x)) * 1 / (1 - e^((-x - log(1 - x)) *  (1/2 (-2 x - x^2 - 2 log(1 - x)))) * (12 x^2 - 12 x^3 - 6 x^4 - 4 x^5 - 3 x^6 - 12 x^2 log(1 - x))/(12 - 12 x^2 + 12 x^3 + 6 x^4 + 4 x^5 + 3 x^6 + 12 x^2 log(1 - x)))
// Tr'(x) = 1/2*(3*x^6 + 4*x^5 + 6*x^4 + 12*x^3 + 12*x^2*log(-x + 1) - 12*x^2)*((3*x^6 + 4*x^5 + 6*x^4 + 12*x^3 + 12*x^2*log(-x + 1) - 12*x^2)*(2*(x + 1/(x - 1) + 1)*(x + log(-x + 1)) + (x^2 + 2*x + 2*log(-x + 1))*(1/(x - 1) + 1))*e^(1/2*(x^2 + 2*x + 2*log(-x + 1))*(x + log(-x + 1)))/(3*x^6 + 4*x^5 + 6*x^4 + 12*x^3 + 12*x^2*log(-x + 1) - 12*x^2 + 12) + 4*(9*x^5 + 10*x^4 + 12*x^3 + 18*x^2 + 12*x*log(-x + 1) - 12*x + 6*x^2/(x - 1))*e^(1/2*(x^2 + 2*x + 2*log(-x + 1))*(x + log(-x + 1)))/(3*x^6 + 4*x^5 + 6*x^4 + 12*x^3 + 12*x^2*log(-x + 1) - 12*x^2 + 12) - 4*(3*x^6 + 4*x^5 + 6*x^4 + 12*x^3 + 12*x^2*log(-x + 1) - 12*x^2)*(9*x^5 + 10*x^4 + 12*x^3 + 18*x^2 + 12*x*log(-x + 1) - 12*x + 6*x^2/(x - 1))*e^(1/2*(x^2 + 2*x + 2*log(-x + 1))*(x + log(-x + 1)))/(3*x^6 + 4*x^5 + 6*x^4 + 12*x^3 + 12*x^2*log(-x + 1) - 12*x^2 + 12)^2)/((3*x^6 + 4*x^5 + 6*x^4 + 12*x^3 + 12*x^2*log(-x + 1) - 12*x^2 + 12)*((3*x^6 + 4*x^5 + 6*x^4 + 12*x^3 + 12*x^2*log(-x + 1) - 12*x^2)*e^(1/2*(x^2 + 2*x + 2*log(-x + 1))*(x + log(-x + 1)))/(3*x^6 + 4*x^5 + 6*x^4 + 12*x^3 + 12*x^2*log(-x + 1) - 12*x^2 + 12) + 1)^2) - 2*(9*x^5 + 10*x^4 + 12*x^3 + 18*x^2 + 12*x*log(-x + 1) - 12*x + 6*x^2/(x - 1))/((3*x^6 + 4*x^5 + 6*x^4 + 12*x^3 + 12*x^2*log(-x + 1) - 12*x^2 + 12)*((3*x^6 + 4*x^5 + 6*x^4 + 12*x^3 + 12*x^2*log(-x + 1) - 12*x^2)*e^(1/2*(x^2 + 2*x + 2*log(-x + 1))*(x + log(-x + 1)))/(3*x^6 + 4*x^5 + 6*x^4 + 12*x^3 + 12*x^2*log(-x + 1) - 12*x^2 + 12) + 1)) + 2*(3*x^6 + 4*x^5 + 6*x^4 + 12*x^3 + 12*x^2*log(-x + 1) - 12*x^2)*(9*x^5 + 10*x^4 + 12*x^3 + 18*x^2 + 12*x*log(-x + 1) - 12*x + 6*x^2/(x - 1))/((3*x^6 + 4*x^5 + 6*x^4 + 12*x^3 + 12*x^2*log(-x + 1) - 12*x^2 + 12)^2*((3*x^6 + 4*x^5 + 6*x^4 + 12*x^3 + 12*x^2*log(-x + 1) - 12*x^2)*e^(1/2*(x^2 + 2*x + 2*log(-x + 1))*(x + log(-x + 1)))/(3*x^6 + 4*x^5 + 6*x^4 + 12*x^3 + 12*x^2*log(-x + 1) - 12*x^2 + 12) + 1))
// Tr''(x) = 1/4*(18144*x^12 + 4032*x^11 + 7200*x^10 + 18432*x^9 + 2304*x^8*(12*log(-x + 1) - 19) - 13824*x^7*(log(-x + 1) - 1) + 6912*x^6*(log(-x + 1) - 9) + 2304*x^5*(12*log(-x + 1) - 5) + 1152*(36*log(-x + 1)^2 - 180*log(-x + 1) + 179)*x^4 - 4608*(18*log(-x + 1)^2 - 54*log(-x + 1) + 37)*x^3 + 6912*(6*log(-x + 1)^2 - 14*log(-x + 1) + 11)*x^2 + 13824*x*(2*log(-x + 1) - 3) - (243*x^24 + 324*x^23*(log(-x + 1) + 4) + 54*(2*log(-x + 1)^2 + 34*log(-x + 1) + 77)*x^22 + 72*(9*log(-x + 1)^2 + 84*log(-x + 1) + 155)*x^21 + 6*(366*log(-x + 1)^2 + 3272*log(-x + 1) + 3401)*x^20 + 4*(2512*log(-x + 1)^2 + 11030*log(-x + 1) + 6895)*x^19 + 4*(324*log(-x + 1)^3 + 7310*log(-x + 1)^2 + 18785*log(-x + 1) + 6963)*x^18 + 16*(378*log(-x + 1)^3 + 3787*log(-x + 1)^2 + 7184*log(-x + 1) + 382)*x^17 + 16*(981*log(-x + 1)^3 + 8010*log(-x + 1)^2 + 6134*log(-x + 1) - 1271)*x^16 + 32*(1602*log(-x + 1)^3 + 5283*log(-x + 1)^2 + 863*log(-x + 1) - 891)*x^15 + 144*(36*log(-x + 1)^4 + 676*log(-x + 1)^3 + 973*log(-x + 1)^2 - 264*log(-x + 1) - 87)*x^14 + 576*(30*log(-x + 1)^4 + 204*log(-x + 1)^3 + 179*log(-x + 1)^2 - 195*log(-x + 1) + 20)*x^13 + 288*(102*log(-x + 1)^4 + 558*log(-x + 1)^3 - 301*log(-x + 1)^2 - 224*log(-x + 1) + 150)*x^12 + 1728*(40*log(-x + 1)^4 + 30*log(-x + 1)^3 - 101*log(-x + 1)^2 + 60*log(-x + 1) + 6)*x^11 + 1728*(4*log(-x + 1)^5 + 38*log(-x + 1)^4 - 63*log(-x + 1)^3 + 10*log(-x + 1)^2 + 88*log(-x + 1) - 59)*x^10 + 6912*(2*log(-x + 1)^5 - log(-x + 1)^4 - 8*log(-x + 1)^3 + 32*log(-x + 1)^2 - 33*log(-x + 1) + 9)*x^9 + 3456*(2*log(-x + 1)^5 - 8*log(-x + 1)^4 + 32*log(-x + 1)^3 - 41*log(-x + 1)^2 + 14*log(-x + 1) + 1)*x^8 + 6912*(2*log(-x + 1)^4 - 2*log(-x + 1)^3 - 3*log(-x + 1)^2 + 4*log(-x + 1) - 1)*x^7 + 6912*(log(-x + 1)^4 - 3*log(-x + 1)^3 + 3*log(-x + 1)^2 - log(-x + 1))*x^6)*e^(x^3 + x^2*log(-x + 1) + 2*x^2 + 4*x*log(-x + 1) + 2*log(-x + 1)^2) + (243*x^24 + 324*x^23*(log(-x + 1) + 4) + 54*(2*log(-x + 1)^2 + 34*log(-x + 1) + 77)*x^22 + 72*(9*log(-x + 1)^2 + 84*log(-x + 1) + 164)*x^21 + 6*(366*log(-x + 1)^2 + 3308*log(-x + 1) + 3743)*x^20 + 4*(2512*log(-x + 1)^2 + 11138*log(-x + 1) + 8137)*x^19 + 4*(324*log(-x + 1)^3 + 7310*log(-x + 1)^2 + 18911*log(-x + 1) + 10380)*x^18 + 16*(378*log(-x + 1)^3 + 3787*log(-x + 1)^2 + 7801*log(-x + 1) + 1175)*x^17 + 16*(981*log(-x + 1)^3 + 8199*log(-x + 1)^2 + 6979*log(-x + 1) - 367)*x^16 + 16*(3204*log(-x + 1)^3 + 10800*log(-x + 1)^2 + 3332*log(-x + 1) + 297)*x^15 + 24*(216*log(-x + 1)^4 + 4056*log(-x + 1)^3 + 5996*log(-x + 1)^2 + 2110*log(-x + 1) - 877)*x^14 + 48*(360*log(-x + 1)^4 + 2448*log(-x + 1)^3 + 3212*log(-x + 1)^2 - 2434*log(-x + 1) + 1695)*x^13 + 48*(612*log(-x + 1)^4 + 3636*log(-x + 1)^3 - 1874*log(-x + 1)^2 + 147*log(-x + 1) + 3592)*x^12 + 192*(360*log(-x + 1)^4 + 294*log(-x + 1)^3 - 771*log(-x + 1)^2 + 1600*log(-x + 1) - 342)*x^11 + 576*(12*log(-x + 1)^5 + 114*log(-x + 1)^4 - 185*log(-x + 1)^3 + 344*log(-x + 1)^2 - 87*log(-x + 1) + 10)*x^10 + 576*(24*log(-x + 1)^5 - 12*log(-x + 1)^4 + 44*log(-x + 1)^3 - 20*log(-x + 1)^2 + 178*log(-x + 1) - 135)*x^9 + 576*(12*log(-x + 1)^5 - 12*log(-x + 1)^4 - 24*log(-x + 1)^3 + 164*log(-x + 1)^2 + 109*log(-x + 1) - 466)*x^8 + 6912*(2*log(-x + 1)^3 + 24*log(-x + 1)^2 - 77*log(-x + 1) + 31)*x^7 + 3456*(20*log(-x + 1)^3 - 79*log(-x + 1)^2 + 20*log(-x + 1) + 20)*x^6 - 6912*(2*log(-x + 1)^3 + 23*log(-x + 1)^2 - 38*log(-x + 1) + 19)*x^5 - 20736*(3*log(-x + 1)^3 - 8*log(-x + 1)^2 + 13*log(-x + 1) - 10)*x^4 - 82944*(log(-x + 1)^2 - 3*log(-x + 1) + 2)*x^3 + 41472*(log(-x + 1)^2 - 2*log(-x + 1) + 1)*x^2)*e^(1/2*x^3 + 1/2*x^2*log(-x + 1) + x^2 + 2*x*log(-x + 1) + log(-x + 1)^2) - 13824*log(-x + 1) + 13824)/(27*x^20 + 54*x^19 + 117*x^18 + 316*x^17 + 2*x^16*(162*log(-x + 1) - 91) + 4*x^15*(54*log(-x + 1) - 11) + 36*x^14*(13*log(-x + 1) - 1) + 72*x^13*(20*log(-x + 1) - 33) + 36*(36*log(-x + 1)^2 - 128*log(-x + 1) + 99)*x^12 - 144*(6*log(-x + 1)^2 - 12*log(-x + 1) - 1)*x^11 + 144*(3*log(-x + 1)^2 + 9*log(-x + 1) - 32)*x^10 + 1728*(log(-x + 1)^2 - 6*log(-x + 1) + 6)*x^9 + 864*(2*log(-x + 1)^3 - 15*log(-x + 1)^2 + 31*log(-x + 1) - 17)*x^8 - 864*(4*log(-x + 1)^3 - 18*log(-x + 1)^2 + 20*log(-x + 1) + 1)*x^7 + 432*(4*log(-x + 1)^3 - 48*log(-x + 1) + 57)*x^6 - 1728*(6*log(-x + 1)^2 - 18*log(-x + 1) + 11)*x^5 + 2592*(2*log(-x + 1)^2 - 2*log(-x + 1) - 3)*x^4 - 5184*x^3*(2*log(-x + 1) - 3) + 1728*x^2*(3*log(-x + 1) - 2) + (27*x^20 + 54*x^19 + 117*x^18 + 316*x^17 + 2*x^16*(162*log(-x + 1) - 91) + 4*x^15*(54*log(-x + 1) - 11) + 36*x^14*(13*log(-x + 1) - 10) + 288*x^13*(5*log(-x + 1) - 9) + 72*(18*log(-x + 1)^2 - 64*log(-x + 1) + 43)*x^12 - 432*(2*log(-x + 1)^2 - 4*log(-x + 1) + 3)*x^11 + 432*(log(-x + 1)^2 - 3*log(-x + 1))*x^10 + 1728*(log(-x + 1)^2 - 5*log(-x + 1) + 5)*x^9 + 864*(2*log(-x + 1)^3 - 15*log(-x + 1)^2 + 30*log(-x + 1) - 17)*x^8 - 1728*(2*log(-x + 1)^3 - 9*log(-x + 1)^2 + 12*log(-x + 1) - 5)*x^7 + 1728*(log(-x + 1)^3 - 3*log(-x + 1)^2 + 3*log(-x + 1) - 1)*x^6)*e^(3/2*x^3 + 3/2*x^2*log(-x + 1) + 3*x^2 + 6*x*log(-x + 1) + 3*log(-x + 1)^2) + 3*(27*x^20 + 54*x^19 + 117*x^18 + 316*x^17 + 2*x^16*(162*log(-x + 1) - 91) + 4*x^15*(54*log(-x + 1) - 11) + 36*x^14*(13*log(-x + 1) - 7) + 360*x^13*(4*log(-x + 1) - 7) + 12*(108*log(-x + 1)^2 - 384*log(-x + 1) + 271)*x^12 - 48*(18*log(-x + 1)^2 - 36*log(-x + 1) + 17)*x^11 + 48*(9*log(-x + 1)^2 - 9*log(-x + 1) - 32)*x^10 + 576*(3*log(-x + 1)^2 - 16*log(-x + 1) + 16)*x^9 + 144*(12*log(-x + 1)^3 - 90*log(-x + 1)^2 + 182*log(-x + 1) - 105)*x^8 - 576*(6*log(-x + 1)^3 - 27*log(-x + 1)^2 + 34*log(-x + 1) - 10)*x^7 + 1728*(log(-x + 1)^3 - 2*log(-x + 1)^2 - 2*log(-x + 1) + 4)*x^6 - 3456*(log(-x + 1)^2 - 3*log(-x + 1) + 2)*x^5 + 1728*(log(-x + 1)^2 - 2*log(-x + 1) + 1)*x^4)*e^(x^3 + x^2*log(-x + 1) + 2*x^2 + 4*x*log(-x + 1) + 2*log(-x + 1)^2) + 3*(27*x^20 + 54*x^19 + 117*x^18 + 316*x^17 + 2*x^16*(162*log(-x + 1) - 91) + 4*x^15*(54*log(-x + 1) - 11) + 36*x^14*(13*log(-x + 1) - 4) + 144*x^13*(10*log(-x + 1) - 17) + 48*(27*log(-x + 1)^2 - 96*log(-x + 1) + 71)*x^12 - 48*(18*log(-x + 1)^2 - 36*log(-x + 1) + 7)*x^11 + 48*(9*log(-x + 1)^2 + 9*log(-x + 1) - 64)*x^10 + 576*(3*log(-x + 1)^2 - 17*log(-x + 1) + 17)*x^9 + 144*(12*log(-x + 1)^3 - 90*log(-x + 1)^2 + 184*log(-x + 1) - 105)*x^8 - 288*(12*log(-x + 1)^3 - 54*log(-x + 1)^2 + 64*log(-x + 1) - 9)*x^7 + 144*(12*log(-x + 1)^3 - 12*log(-x + 1)^2 - 84*log(-x + 1) + 109)*x^6 - 576*(12*log(-x + 1)^2 - 36*log(-x + 1) + 23)*x^5 + 864*(4*log(-x + 1)^2 - 6*log(-x + 1) - 1)*x^4 - 1728*x^3*(2*log(-x + 1) - 3) + 1728*x^2*(log(-x + 1) - 1))*e^(1/2*x^3 + 1/2*x^2*log(-x + 1) + x^2 + 2*x*log(-x + 1) + log(-x + 1)^2) - 3456*x + 1728)

 // derivative((12*x^2 - 12*x^3-6*x^4-4*x^5-3*x^6-12*x^2*log(1-x))/(12-12*x^2+12*x^3+6*x^4+4*x^5+3*x^6+12*x^2*log(1-x))*1/(1-e^((-x-log(1-x))*(1/2*(-2*x-x^2 -2*log(1-x))))*(12*x^2-12*x^3-6*x^4-4*x^5-3*x^6-12*x^2*log(1-x))/(12-12*x^2+12*x^3+6*x^4+4*x^5+3*x^6+12*x^2*log(1-x))) ,x)
