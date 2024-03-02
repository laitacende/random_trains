const express = require('express');
const path = require('path');
const app = express();
const port = 3000;

const modified = require("./boltzmannSamplers.js");
const orginal = require("./boltzmannSamplers_orginal");

app.use(express.static(__dirname));
app.set('view engine', 'ejs');

app.get('/', (req, res) => {
    res.render('index', {
        // msg: JSON.stringify(tr)
        //     msg: JSON.stringify(modified.gammaTrain(0.3))
       msg: JSON.stringify(orginal.gammaTrain(0.25))
    });
});

app.listen(3000, () => console.log('Listening on port 3000.'));


var tr = {
    "wagons": [
        {
            "wagon": [
                {
                    "planks": [
                        "○",
                        "○"
                    ],
                    "wheels": [
                        [
                            "•",
                            "•"
                        ]
                    ]
                },
                {
                    "planks": [
                        "○",
                        "○"
                    ],
                    "wheels": [
                        [
                            "•",
                            "•"
                        ]
                    ]
                }
            ],
            "passengers": [{"passenger": [["p", "p", "p"], ["p", "p", "p"]]},
                {"passenger": [["p"], ["p", "p"]]}]
        },
        {
            "wagon": [
                {
                    "planks": [
                        "○",
                        "○"
                    ],
                    "wheels": [
                        [
                            "•",
                            "•"
                        ]
                    ]
                },
            ],
            "passengers": []
        },
        {
            "wagon": [
                {
                    "planks": [
                        "○",
                        "○"
                    ],
                    "wheels": [
                        [
                            "•",
                            "•"
                        ]
                    ]
                },
            ],
            "passengers": []
        }
    ]
};