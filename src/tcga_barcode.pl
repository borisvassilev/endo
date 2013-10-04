/*
Copyright (c) 2013 Boris Vassilev, University of Helsinki

This file is part of "endo" Project. It is subject to the license
terms in the LICENSE file found in the top-level directory of
this distribution and at http://opensource.org/licenses/MIT. No
part of "endo" Project, including this file, may be copied,
modified, propagated, or distributed except according to the
terms contained in the LICENSE file.
*/

% This module provides a DCG for parsing TCGA barcodes and a
% predicate for extracting information from them. Some validation
% of the contents of the barcode takes place.
:- module(tcga_barcode,
    [   barcode//2,
        extract_barcode/4
    ]).

:- use_module(library('dcg/basics')).

%! extract_barcode(ExtractType, Type, Content, Extracted) is det
%
% Extract the information relevant for a barcode of type
% _ExtractType_ from the barcode of type _Type_ with content
% _Content_. Result is a compound term named after the extracted
% information.
extract_barcode(participant, Type, [TSS,P|_], participant(TSS,P)) :-
    Type \== tss.

% Barcode DCG
% -----------
%
% This DCG attempts to parse a string representing a valid TCGA
% barcode as defined here:
% https://wiki.nci.nih.gov/display/TCGA/Working+with+TCGA+Data
barcode(Type, Contents) -->
    "TCGA-",
    tss(Type, Contents),
    !. % make deterministic

tss(Type, [TSS|More]) -->
    string_without("-", TSSChars),
    {   % validation
        atom_codes(TSS, TSSChars),
        tcga_data:tissue_source_site(TSS, _, _, _)
    },
    tss_(Type, More).
tss_(Type, More) -->
    "-", 
    participant(Type, More).
tss_(tss, []) --> [].

participant(Type, [P|More]) -->
    alphanumeric_string(ParticipantChars),
    {   atom_chars(P, ParticipantChars)
    },
    participant_(Type, More).
participant_(Type, More)  -->
    [0'-,First],
    what_sample(First, Type, More).
participant_(participant, []) --> [].

% What kind of sample is this?
what_sample(0'C, drug, ['C'-Drug]) -->
    drug(Drug).
what_sample(0'D, drug, ['D'-Drug]) -->
    drug(Drug).
what_sample(0'H, drug, ['H'-Drug]) -->
    drug(Drug).
what_sample(0'I, drug, ['I'-Drug]) -->
    drug(Drug).
what_sample(0'T, drug, ['T'-Drug]) -->
    drug(Drug).
what_sample(0'E, examination, [Examination]) -->
    examination(Examination).
what_sample(0'S, surgery, [Surgery]) -->
    surgery(Surgery).
what_sample(0'R, radiation, [Radiation]) -->
    radiation(Radiation).
% sample
what_sample(First, Type, More) -->
    {   % validaton
        code_type(First, digit)
    },
    sample(First, Type, More).

drug(Drug) -->
    integer(Drug).
examination(Examination) -->
    integer(Examination).
surgery(Surgery) -->
    integer(Surgery).
radiation(Radiation) -->
    integer(Radiation).

sample(First, Type, [Sample|More]) -->
    digits(Rest),
    {   % validation
        number_codes(Sample, [First|Rest]),
        tcga_data:sample_type(Sample, _, _)
    },
    sample_(Type, More).
sample_(Type, More) -->
    [Vial],
    {   code_type(Vial, upper)
    }, !,
    portion(Vial, Type, More).
sample_(sample, []) --> [].

portion(Vial, Type, [V-N|More]) -->
    "-", digit(D1), digit(D2),
    {   char_code(V, Vial),
        atom_codes(N, [D1,D2])
    },
    portion_(Type, More).
portion_(Type, More) -->
    "-",
    slide(Type, More).
portion_(Type, [Analyte|More]) -->
    [AnalyteChar],
    {   %validation
        char_code(Analyte, AnalyteChar),
        tcga_data:portion_analyte(Analyte, _)
    },
    analyte(Type, More).
portion_(portion, []) --> [].

slide(slide, [TSID-Slide]) -->
    [S1,S2,S],
    {   % validation
        atom_codes(TSID, [S1,S2]),
        slide_id(TSID),
        code_type(S, alnum),
        char_code(Slide, S)
    }.

slide_id('TS'). % top
slide_id('BS'). % bottom
slide_id('MS'). % middle slide

analyte(Type, More) -->
    aliquot(Type, More).
analyte(analyte, []) --> [].

aliquot(aliquot, [Aliquot,Center]) -->
    [0'-,A1,A2,A3,A4,0'-],
    {   % validation
        maplist(type_code(alnum), [A1,A2,A3,A4]),
        atom_codes(Aliquot, [A1,A2,A3,A4])
    },
    center(Center).
% NOTE: also used by not_code/2 predicate from this file
type_code(Type, Code) :-
    code_type(Code, Type).

center(Center) -->
    nonblanks(CenterChars),
    {   % validation
        atom_chars(Center, CenterChars),
        tcga_data:center_code(Center, _, _, _, _)
    }.

% Additional dcg basics
% ---------------------
% Read a string of upper case letters deterministically.
upper_string([U|T]) -->
    [U],
    {   code_type(U, upper)
    }, !,
    upper_string(T).
upper_string([]) --> [].

% Read a string of alphanumeric characters deterministically
alphanumeric_string([C|T]) -->
    [C],
    {   code_type(C, alnum)
    }, !,
    alphanumeric_string(T).
alphanumeric_string([]) --> [].

% UNIT TESTS
% ----------
%
% Test compliance of the barcode DCG to the somewhat informal
% specification provided at
% https://wiki.nci.nih.gov/display/TCGA/Working+with+TCGA+Data
%
% As values defined in the Code Tables provided by the TCGA are
% used to validate the barcode while parsing it, those are
% extensively tested.
%
% The test names are as descriptive as possible and further
% documentation should not be necessary.
%
% Helper predicates for unit tests
% --------------------------------
%
tissue_source_site_outside_domain([F,S]) :-
    code([digit,upper], F),
    code([digit,upper], S),
    atom_codes(TSS, [F,S]),
    \+ tcga_data:tissue_source_site(TSS, _, _, _).

participant_outside_domain([C1,C2,C3,C4]) :-
    select(C, [C1,C2,C3,C4], Rest),
    Rest = [0'Z,0'1,0'B],
    not_code([alnum], C).

portion_domain(P) :-
    code([digit], F),
    code([digit], S),
    atom_codes(P, [F,S]).

slide_domain(ID-S) :-
    member(ID, ['TS', 'BS', 'MS']),
    code([alnum], C),
    char_code(S, C).

portion_analyte_outside_domain([C]) :-
    code([upper], C),
    char_code(PA, C),
    \+ tcga_data:portion_analyte(PA, _).

plate_outside_domain(P) :-
    participant_outside_domain(P).

center_code_outside_domain([F,S]) :-
    code([digit], F),
    code([digit], S),
    atom_codes(C, [F,S]),
    \+ tcga_data:center_code(C, _, _, _, _).

% Enumerate all codes that belong to a list of (non-overlapping)
% ranges in the form FromCode-ToCode or a list of types
% recognized by code_type/2 SWI-Prolog built-in
code(Specs, Code) :-
    member(S, Specs), % for each specification
    code_(S, Code). % enumerate fitting codes
% A range
code_(From-To, Code) :- !,
    between(From, To, Code).
% A type
code_(Type, Code) :- !,
    code_type(Code, Type).

% Enumerate all ASCII codes between 1 and 127 that _do not_
% belong to a list of ranges FromCode-ToCode or a type recognized
% by code_type/2
not_code(Specs, C) :-
    numlist(0, 127, Codes), % all possible codes
    % use the specifications to progressively filter
    foldl( not_code_, Specs, Codes, NotCodes ),
    member(C, NotCodes). % enumerate all codes left
not_code_(From-To, Codes, NotCodes) :- !,
    exclude( between(From, To), Codes, NotCodes ).
not_code_(Type, Codes, NotCodes) :- !,
    exclude( type_code(Type), Codes, NotCodes ). % notice type_code/2

% Tests for TCGA-Barcode DCG
% --------------------------
%
:- begin_tests(barcode_tests, [setup(tcga_data:load_data)]).

test(tcga_only_code_fail, 
    [   fail
    ]
) :-
    phrase( barcode(_T, _C), "TCGA" ).
    
test(tcga_not_tcga_fail, 
    [   fail
    ]
) :-
    phrase( barcode(_T, _C), "TSGA-02" ).
    
test(tcga_not_tss_fail, 
    [   fail
    ]
) :-
    phrase( barcode(_T, _C), "TSGA-NOTTSS" ).
    
test(tss_all,
    [   forall(tcga_data:tissue_source_site(TSS, _, _, _)),
        true(T-C == tss-[TSS])
    ]
) :-
    atom_codes(TSS, TSSCodes),
    append("TCGA-", TSSCodes, Barcode),
    phrase( barcode(T, C), Barcode ).

test(tss_not_tss_fail,
    [   fail,
        forall(tissue_source_site_outside_domain(TSS))
    ]
) :-
    append("TCGA-", TSS, Barcode),
    phrase( barcode(_T, _C), Barcode ).

test(participant,
    [   true(T-C == participant-['02','0001'])
    ]
) :-
    phrase( barcode(T, C), "TCGA-02-0001" ).

test(participant_not_alphanumeric_fail,
    [   fail,
        forall(participant_outside_domain(P))
    ]
) :-
    append("TCGA-02-", P, Barcode),
    phrase( barcode(participant, _C), Barcode ).

test(drugc,
    [   true(T-C == drug-['02','0001','C'-1])
    ]
) :-
    phrase( barcode(T, C), "TCGA-02-0001-C1" ).

test(drugd,
    [   true(T-C == drug-['02','0001','D'-2])
    ]
) :-
    phrase( barcode(T, C), "TCGA-02-0001-D2" ).

test(drugh,
    [   true(T-C == drug-['02','0001','H'-33])
    ]
) :-
    phrase( barcode(T, C), "TCGA-02-0001-H33" ).

test(drugi,
    [   true(T-C == drug-['02','0001','I'-666])
    ]
) :-
    phrase( barcode(T, C), "TCGA-02-0001-I0666" ).

test(drugt,
    [   true(T-C == drug-['02','0001','T'-999999])
    ]
) :-
    phrase( barcode(T, C), "TCGA-02-0001-T999999" ).

test(examination,
    [   true(T-C == examination-['02','0001',3124])
    ]
) :-
    phrase( barcode(T, C), "TCGA-02-0001-E3124" ).

test(surgery,
    [   true(T-C == surgery-['02','0001',145])
    ]
) :-
    phrase( barcode(T, C), "TCGA-02-0001-S145" ).

test(radiation,
    [   true(T-C == radiation-['02','0001',2])
    ]
) :-
    phrase( barcode(T, C), "TCGA-02-0001-R2" ).

test(sample_all,
    [   forall(tcga_data:sample_type(S, _, _)),
        true(sample-C == sample-['02','0001',S])
    ]
) :-
    number_codes(S, Codes),
    append("TCGA-02-0001-", Codes, Barcode),
    phrase( barcode(sample, C), Barcode ).

test(portion_all,
    [   forall(portion_domain(P)),
        true(portion-C == portion-['02','0001',1,'C'-P])
    ]
) :-
    atom_codes(P, PCodes),
    append("TCGA-02-0001-01C-", PCodes, Barcode),
    phrase( barcode(portion, C), Barcode ).

test(slide,
    [   forall( slide_domain(ID-S) ),
        true(T-C == slide-['02','0001',1,'C'-'01', ID-S])
    ]
) :-
    atom_codes(ID, IDCodes),
    atom_codes(S, SCodes),
    append(["TCGA-02-0001-01C-01-", IDCodes, SCodes], Barcode),
    phrase( barcode(T, C), Barcode ).


test(slide_id_fail,
    [   fail 
    ]
) :-
    phrase( barcode(_T, _C), "TCGA-02-0001-01C-01-AS1" ).

test(slide_fail,
    [   fail 
    ]
) :-
    phrase( barcode(_T, _C), "TCGA-02-0001-01C-01-BSAA" ).

test(analyte_all,
    [   forall( tcga_data:portion_analyte(PA, _) ),
        true(T-C == analyte-['02','0001',1,'C'-'01',PA])
    ]
) :-
    atom_codes(PA, PACodes),
    append("TCGA-02-0001-01C-01", PACodes, Barcode),
    phrase( barcode(T, C), Barcode ).

test(analyte_fail,
    [   forall( portion_analyte_outside_domain(PA) ),
        fail
    ]
) :-
    append("TCGA-02-0001-01C-01", PA, Barcode),
    phrase( barcode(analyte, _C), Barcode ).

test(plate_fail,
    [   forall(plate_outside_domain(P)),
        fail
    ]
) :-
    append(["TCGA-02-0001-01C-01D-", P, "-01"], Barcode),
    phrase( barcode(_T, _C), Barcode ).

test(center_all,
    [   forall( tcga_data:center_code(Center, _, _, _, _) ),
        true(T-C == aliquot-['02','0001',1,'C'-'01','D','0182',Center])
    ]
) :-
    atom_codes(Center, CenterCodes),
    append("TCGA-02-0001-01C-01D-0182-", CenterCodes, Barcode),
    phrase( barcode(T, C), Barcode ).

test(center_fail,
    [   forall( center_code_outside_domain(Center) ),
        fail
    ]
) :-
    append("TCGA-02-0001-01C-01D-0182-", Center, Barcode),
    phrase( barcode(_T, _C), Barcode ).

:- end_tests(barcode_tests).

