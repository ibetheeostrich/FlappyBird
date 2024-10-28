function usydcolours()
%% define colours
% some of the default matlab colours are quite bad so I decided to overwrite them
% the new scheme on the recent matlab versions is a lot better but I still prefer to start
% off with black, then blue or red, ...
% especially if you only have one set of data it can save you a lot in colour printing
AccentBlue	            =    [1	    72	164	]/255;
AccentGrey	            =    [241	241	241	]/255;
AccentYellow	        =    [255	184	0	]/255;
MasterbrandBlack	    =    [10	10	10	]/255;
MasterbrandCharcoal	    =    [66	66	66	]/255;
MasterbrandOchre	    =    [230	70	38	]/255;
MasterbrandWhite	    =    [255	255	255	]/255;
SecondaryBeige	        =    [253	202	144	]/255;
SecondaryBlue	        =    [78	152	211	]/255;
SecondaryDarkGreen	    =    [0	    126	59	]/255;
SecondaryDarkSeafoam	=    [0	    164	133	]/255;
SecondaryIvory	        =    [248	239	221	]/255;
SecondaryLemon	        =    [251	243	141	]/255;
SecondaryLightBlue	    =    [145	189	229	]/255;
SecondaryLightGreen	    =    [189	220	150	]/255;
SecondaryLightPink	    =    [248	185	204	]/255;
SecondaryLightSeafoam	=    [104	198	182	]/255;
SecondaryLilac	        =    [184	150	198	]/255;
SecondaryMaroon	        =    [122	32	0	]/255;
SecondaryOrange	        =    [249	161	52	]/255;
SecondaryPeach	        =    [247	156	114	]/255;
SecondaryPink	        =    [214	81	157	]/255;
SecondaryPurple	        =    [127	63	152	]/255;
DarkEucalypt	        =    [37	88	77	]/255;
Eucalypt	            =    [113	164	153	]/255;
NeutralGrey	            =    [224	224	224	]/255;
LightOchre	            =    [255	173	140	]/255;
UpdatedOchre	        =    [231	71	38	]/255;

% if you want to keep using the standard matlab color remove this
set(groot,'defaultAxesColorOrder',[UpdatedOchre; AccentBlue; AccentYellow; SecondaryDarkGreen]);
end

