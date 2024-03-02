var canvas = document.getElementById("canvas");
var canvasContext = canvas.getContext("2d");

canvas.width = 0.9 * window.innerWidth;     // equals window dimension
canvas.height = 0.9 * window.innerHeight;
console.log(msgTr);
function drawCircle(color, radius, x, y) {
    canvasContext.beginPath();
    canvasContext.arc(x, y, radius, 0, 2 * Math.PI, false);
    canvasContext.fillStyle = color;
    canvasContext.fill();
}

function drawEllipse(color, radius, x, y) {
    canvasContext.beginPath();
    canvasContext.ellipse(x, y, radius * 0.5, radius, 0, 0, Math.PI*2, true);
    canvasContext.fillStyle = color;
    canvasContext.fill();
}

// jaki procent canvasa może zająć szerokość
// x, y pozycje gdzie zaczac rysowac
function drawWagon(wagon, width, x, y) {
   // let width = canvas.width * percent;
    canvasContext.fillStyle = "#79b2f2";
    canvasContext.fillRect(x, y, width, width * 0.5);
    // rozdziel na planki z kołami i bez
    let withWheel = [];
    let withoutWheel = []
    for (let i = 0; i < wagon.length; i++) {
        if (wagon[i]["wheels"].length != 0) {
            withWheel.push(wagon[i]);
        } else {
            withoutWheel.push(wagon[i]);
        }
    }
    // narysuj planki bez kół do góry
    let step = width / withoutWheel.length;
    let currX = x;
    for (let i = 0; i < 2 * withoutWheel.length; i++) {
        drawCircle("#234b78", width * 0.02, currX, y + width * 0.02 / 2);
        currX += step / 2;
    }

    // z kołami narysuj na dole
    step = width / withWheel.length;
    currX = x;
    for (let i = 0; i < 2 * withWheel.length; i++) {
        drawCircle("#234b78", width * 0.02, currX, y + width * 0.5);
        currX += step / 2;
    }

    // narysuj koła na środku każdego z planków
    // każdy ma jedno koło, ale mogą różnić się liczbą punktów
    currX = x + step / 4;
    for (let i = 0; i < withWheel.length; i++) {
        // baza koła
        drawCircle("#858a8f", width * 0.07, currX, y + width * 0.57);
        canvasContext.fillStyle = "black";
        canvasContext.font = "1em Arial";
        canvasContext.fillText(withWheel[i].wheels[0].length, currX -  width * 0.07 / 2 , y + width * 0.57);
        currX += step;
    }

}

function drawPassenger(pas, radius, x, y) {
    // głowa jest pierwsza
    drawCircle("#eddf45", radius * 0.3, x, y);
    canvasContext.fillStyle = "black";
    canvasContext.font = "1em Arial";
    canvasContext.fillText(pas[0].length, x , y + radius * 0.2);

    // tłów
    drawEllipse("#eddf45", radius*0.7, x, y + radius*1);
    canvasContext.fillStyle = "black";
    canvasContext.font = "1em Arial";
    canvasContext.fillText(pas[1].length, x , y + 1*radius);

}

function drawTrain(train) {
    // narysuj wagony
    let noWagons = train["wagons"].length;
    let smallerCanvW = canvas.width - 60;
    let widthWagon = smallerCanvW / noWagons * 0.8;
    let widthSpace = smallerCanvW * 0.2 / (noWagons - 1);
    let currX = 30;
    for (let i = 0; i < noWagons; i++) {
        drawWagon(train["wagons"][i]["wagon"], widthWagon, currX, 0.4 * canvas.height);
        // draw passengers
        // TODO USUN TO
        if (train["wagons"][i]["passengers"].length != 0) {
            let noPass = train["wagons"][i]["passengers"].length;
            let widthPass = widthWagon / noPass * 0.9 * 0.2;
            let widthSpacePass = 0;
            let currXPass = currX;
            let rad = widthPass;
            for (let j = 0; j < noPass; j++) {
                drawPassenger(train["wagons"][i]["passengers"][j]["passenger"], rad,
                    currXPass + rad, 0.45 * canvas.height);
                currXPass += widthPass + widthSpacePass;
            }
        }
        currX +=  widthWagon + widthSpace;
    }
}
drawTrain(JSON.parse(msgTr));
// drawTrain({
// "wagons": [
//      {
//          "wagon": [
//               {
//                  "planks": [
//                          "○",
//                          "○"
//                          ],
//                          "wheels": [
//                              [
//                                  "•",
//                                  "•"
//                              ]
//                          ]
//                 },
//                  {
//                      "planks": [
//                              "○",
//                              "○"
//                          ],
//                          "wheels": [
//                                   [
//                                      "•",
//                                      "•"
//                                     ]
//                           ]
//                  }
//              ],
//              "passengers": [{"passenger": [["p", "p", "p"], ["p", "p", "p"]]},
//                  {"passenger": [["p"], ["p", "p"]]}]
//          },
//     {
//         "wagon": [
//             {
//                 "planks": [
//                     "○",
//                     "○"
//                 ],
//                 "wheels": [
//                     [
//                         "•",
//                         "•"
//                     ]
//                 ]
//             },
//         ],
//         "passengers": []
//     },
//     {
//         "wagon": [
//             {
//                 "planks": [
//                     "○",
//                     "○"
//                 ],
//                 "wheels": [
//                     [
//                         "•",
//                         "•"
//                     ]
//                 ]
//             },
//         ],
//         "passengers": []
//     }
//   ]
// });

// drawWagon([{"planks": [
//             "○",
//             "○"
//         ],
//         "wheels": [
//             [
//                 "•",
//                 "•"
//             ]
//         ]},
//     {"planks": [
//             "○",
//             "○"
//         ],
//         "wheels": [
//         ]},
//     {"planks": [
//             "○",
//             "○"
//         ],
//         "wheels": [
//             [
//                 "•",
//                 "•",
//                 "•",
//                 "•"
//             ]
//         ]}], 0.2, 50, 50);

// drawWagon([{"planks": [
//         "○",
//         "○"
//     ],
//     "wheels": [
//     ]}], 0.2, 50, 50);