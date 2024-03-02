const MersenneTwister = require('mersenne-twister');
var betainc = require( '@stdlib/math-base-special-betainc' );
let generator = new MersenneTwister();
// generator.random_incl();  on [0, 1]
let loopsBound = -1;
function uniform() {
    return generator.random_incl();
}

function bernoulli(p) {
    let rand = uniform();
    if (rand <= p) {
        return 1;
    } else {
        return 0;
    }
}

// TODO odrzuc jesli wieksze niz end
function geometric(p, start, end) {
    let loops = 0;
    let rand = uniform();
    let tmp = 0.0;
    let k = start;
    while (rand >= tmp) {
        tmp += Math.pow(p, k) * (1 - p);
        k += 1;
        if (loopsBound > 0 && loops >= loopsBound) {
            return -1;
        }
        loops += 1;
    }
    return k;
}

function geometricTruncated(p, start) {
    let loops = 0;
    let rand = uniform();
    let tmp = 0.0;
    let k = start;
    let tmpPrev = 1;
    let cdfScale = 1 - Math.pow(p, Math.floor(start));
    cdfScale = 1 - cdfScale;
    while (rand >= tmp && Math.abs(tmp- tmpPrev) > 0.0000001) {
        // console.log(k + " " + tmp + " " + tmpPrev + " " + rand + " " + Math.abs(tmp- tmpPrev))
        tmpPrev = tmp;
        tmp += Math.pow(p, k)  / cdfScale * (1-p);
        k += 1;
        if (loopsBound > 0 &&  loops >= loopsBound) {
            return -1;
        }
        loops += 1;
    }
    return k - 1;
}

function loga(p, start) {
    let loops = 0;
    let rand = uniform();
    let tmp = 0.0;
    let k = start;
    while (rand >= tmp) {
        tmp += Math.pow(p, k) / k * -1.0 / Math.log(1 - p);
        k += 1;
        if (loopsBound > 0 && loops >= loopsBound) {
            return -1;
        }
        loops += 1;
    }
    return k;
}


function logaTruncated(p, start) {
    let loops = 0;
    let rand = uniform();
    let tmp = 0.0;
    let cdfScale = 1.0 + betainc(p, start, 0.0000000000000001, false) / Math.log(1.0 - p);
    cdfScale = 1 - cdfScale;
    let k = start;
    while (rand >= tmp) {
        tmp += Math.pow(p, k) / k * -1.0 / Math.log(1 - p) / cdfScale;
        k += 1;
        if (loopsBound > 0 && loops >= loopsBound) {
            return -1;
        }
        loops += 1;
    }
    return k-1;
}

function factorial(n) {
    let p = BigInt(1)
    for (let i = BigInt(n); i > 0; i--) p *= i
    return p
}

function poisson(p) {
    let rand = uniform();
    let tmp = 0.0;
    let k = 0;
    let fact = 1;
    while (rand >= tmp) {
        tmp += Math.pow(p, k) * Math.exp(-p) / fact;
        k += 1;
        fact *= k;
    }
    return k;
}

module.exports = {geometric, loga, bernoulli, poisson, logaTruncated, geometricTruncated};