(* ::Package:: *)

(* --- HOW TO IMPORT THIS PACKAGE INTO A NOTEBOOK --- *)
(* Get["<path>/dmbog_std_mathematica.wl"];            *)

(* --- TIPS ---                                                                            *)
(* > Evaluate notebook to get Mathematica to suggest autocompletion for `Std\[FilledSquare]...` funtions *)
(* > To type `\[FilledSquare]` symbol: [press Esc] -> [type `fsq`] -> [press Esc]                          *)

(* --- LAST UPDATED --- *)
(* 2024.08.24           *)


(* Standard styles for different level headers *)
Std\[FilledSquare]Header1[text_] := Print@Framed@Style["--- "<>text<>" ---", Bold, FontSize -> 26];
Std\[FilledSquare]Header2[text_] := Print@Framed@Style[ "-- "<>text<>" --",  Bold, FontSize -> 22];
Std\[FilledSquare]Header3[text_] := Print@Framed@Style[  "- "<>text<>" -",   Bold, FontSize -> 18];
Std\[FilledSquare]Header4[text_] := Print@Framed@Style[   " "<>text<>" ",    Bold, FontSize -> 14];

(* Standart parser for solutions saved as JSON series                                                             *)
(* Takes directory with `params.json` (optional) and iterations saved as JSON and enumerated with leading zeroes. *)
Std\[FilledSquare]ParseSolutionAsJSONSeries[dir_] := Module[{filenamesParams, filenamesIterations, params, iterations},
	Std\[FilledSquare]Header1["Parsing solution as 'JSON Series'"];
	
	(* Count files in solution directory *)
	filenamesParams     = FileNames["params.json",                    dir];
	filenamesIterations = FileNames[RegularExpression["\\d+\\.json"], dir];
	
	Print@Row@{"Expected files -> {0000.json, ... , 0XXX.json, params.json}"};
	Print@Row@{"Found parameter files -> ", FileNameTake/@filenamesParams};
	Print@Row@{"Found iteration files -> ", FileNameTake/@filenamesIterations};
	
	(* Parse problem params (if present) *)
	params = If[Length@filenamesParams > 0, Import[filenamesParams[[1]]], {}];
	Print@Row@{"Solution params -> ", Grid[{{"--- KEY ---"~Style~Bold, "--- VALUE ---"~Style~Bold}}~Join~Transpose@{Keys@params, Values@params}, Frame -> All]};
	
	(* Parse problem iterations (if present) *)
	iterations = If[Length@filenamesIterations > 0, Table[Import[filenamesIterations[[k]]], {k, 1, Length@filenamesIterations}], {}];
	Print@Row@{"Saved iteration count -> ", Length@filenamesIterations};
	
	Return[{params, iterations}];
];


(* Standard style preset for 1D plots *)
Std\[FilledSquare]PlotStyle = {
	PlotStyle      -> {Directive[Red, Thick], Directive[Blue, Thick], Directive[Purple, Thick]},
	AxesStyle      -> Directive[Black, Thick],
	LabelStyle     -> Directive[Black, Bold, FontSize -> 14],
	TicksStyle     -> Directive[Black, Bold, FontSize -> 10],
	GridLines      -> {Automatic, None},
	GridLinesStyle -> Directive[Gray, Dashed],
	PlotRange      -> Full,
	ImageSize      -> 500
};


(* "Smart" relative error function that appropriately handles annoying cases: where solution is = 0 and cases where error is precisely 0 *)
(* CASE: approx == precise            => relative error = 0                                                                               *)
(* CASE: approx == 0 || precise == 0   => relative error uses their average to normalize                                                   *)
Std\[FilledSquare]RelativeError[approx_, precise_] := If[Abs@approx + Abs@precise != 0, (approx-precise)/((Abs@approx+Abs@precise)/2), 0];

(* Same thing but for lists: computes error for each pair {approx\[LeftDoubleBracket]k\[RightDoubleBracket], precise\[LeftDoubleBracket]k\[RightDoubleBracket]} and returns result as a list *)
Std\[FilledSquare]RelativeListError[approx_List, precise_List] := Table[Std\[FilledSquare]RelativeError[approx[[elem]], precise[[elem]]], {elem, 1, Length@approx}];


(* Frames given expression, if list is passed all arguments get positioned in a centered column *)
(* Mostly useful when drawing PlotLabel with additional info (especially inside Manipulate[])   *)
Std\[FilledSquare]FramedLabel[lines_List] := Framed[Column[lines, Center]];
Std\[FilledSquare]FramedLabel[line_]      := Std\[FilledSquare]FramedLabel[{line}];


Print@Style["Imported dmbog_std_mathematica.wl", Plain, FontSize -> 16, Darker@Purple];
