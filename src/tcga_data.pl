/*
Copyright (c) 2013 Boris Vassilev, University of Helsinki

This file is part of "endo" Project. It is subject to the license
terms in the LICENSE file found in the top-level directory of
this distribution and at http://opensource.org/licenses/MIT. No
part of "endo" Project, including this file, may be copied,
modified, propagated, or distributed except according to the
terms contained in the LICENSE file.
*/

% This module provides predicates for loading all relevant data
% to the Prolog database. Some validation of the input and the
% correctness of the parsing process is provided.
:- module(tcga_data,
    [load_data/0
    ]).

:- use_module(library(csv)).
:- use_module(tcga_barcode).

% TCGA Code Tables
% ----------------
%
% This is a list of the code table files and the expected columns
% in each of these files. The functor name is used for the name
% of the fact to be added to the database, while the arguments
% must match the type of the converted values within each row.
tcga_code_tables(
    [   code_table(
            '../data/tissueSourceSite.tsv',
            % TSS Code, Source Site, Study Name, BCR
            tissue_source_site(atom, atom, atom, atom)
        ),
       code_table(
            '../data/sampleType.tsv',
            % Code, Definition, Short Letter Code
            sample_type(integer, atom, atom)
        ),
        code_table(
            '../data/portionAnalyte.tsv',
            % Code, Definition
            portion_analyte(atom, atom)
        ),
        code_table(
            '../data/centerCode.tsv',
            % Code, Center Name, Center Type, Display Name, Short Name
            center_code(atom, atom, atom, atom, atom)
        )
    ]
).

% TCGA BRCA gene expression file
% ------------------------------
%
% This is the file that contains the normalized gene expression
% data. The assumed structure of this file is a header row with
% the barcodes of the samples, followed by one row for each gene
% expression profile for all patients.
tcga_brca_ge_file(
    '../data/brca_ge_norm.tsv'
    %'../data/foo.tsv' % a dummy file smaller than the working file
).

load_data :-
    format('loading TCGA code tables '), flush_output,
    time( load_tcga_code_tables ),
    format(' done~n'), flush_output,
    format('loading BRCA gene expression data '), flush_output,
    time( load_tcga_brca_ge_file ),
    format(' done~n'), flush_output.

% Code Tables read and load
% -------------------------
%
% Read the Code Tables provided by the TCGA and load them as
% facts to the database. These are used to get more information
% about the samples. It allows some level of validation for the
% barcodes of all analyzed samples.
%
%! load_tcga_code_tables is det
% Loads all code tables listed in `tcga_code_tables(CodeTable)` to
% the database.
load_tcga_code_tables :-
    tcga_code_tables(CodeTables),
    forall(
        member(code_table(File, Functor), CodeTables),
        (   functor(Functor, Name, Arity),
            abolish(Name/Arity),
            dynamic(Name/Arity),
            load_code_table_file(File, Functor, Name, Arity, 1),
            format('.'), flush_output
        )
    ).

%! load_code_table_file(+File, +Functor, +Name, +Arity, +Skip) is det
% Loads one of the code tables to the database, making use of the
% standard library(csv). 
load_code_table_file(File, Functor, Name, Arity, Skip) :-
    FirstLine is Skip + 1,
    forall(
        csv_read_file_row(
            File,
            R,
            [line(Line),arity(Arity),separator(0'\t),convert(false)]
        ),
        (   Line >= FirstLine
        ->  validate_table_row(R, Functor, Name, Arity, Row),
            assertz(Row)
        ;   true
        )
    ).

%! validate_table_row(+R, +Functor, +Name, +Arity, -Row) is det
% Validate the arguments of the functor returned by library(csv)
% and create a functor with converted values
validate_table_row(R, Functor, Name, Arity, Row) :-
    functor(Row, Name, Arity),
    foreach(
        (   arg(Arg, Functor, Type),
            arg(Arg, R, RawValue),
            convert_field(Type, RawValue, Value)
        ),
        arg(Arg, Row, Value)
    ).

%! convert_field(+Type, +Field, -Converted) is det
% Convert one field to the given type. Also validates the input.
convert_field(atom, Field, Field).
convert_field(number, Field, Number) :-
    atom_number(Field, Number).
convert_field(integer, Field, Int) :-
    atom_number(Field, Int),
    integer(Int).
convert_field(float, Field, Float) :-
    atom_number(Field, Float),
    float(Float).

% Gene Expression Data read and load
% ----------------------------------
%
% The gene expression data contains all normalized expression
% values for over 500 patients and 18 thousand genes. This is a
% big file and loading it takes a while (more than a minute on a
% relatively low-power laptop).
%
%! load_tcga_brca_ge_file is det
% Load the breast cancer patients gene expression file.
load_tcga_brca_ge_file :-
    tcga_brca_ge_file(FileName),
    setup_call_cleanup(
        open(FileName, read, File),
        load_ge_file(File),
        close(File)
    ).

%! load_ge_file(+File) is det
% Read the header line separately. This is necessary because this
% line contains the sample names. All following lines contain
% gene expression profiles.
load_ge_file(File) :-
    read_header_line(File),
    read_expression_profile_lines(File).

%! read_header_line(+File) is det
% Read the first line of the file, assuming it contains the
% barcodes for all samples. Parse the barcode and extract the
% participant. Add facts to the database with the name of the
% sample as provided, the extracted information about the
% participant and the position of the sample within the row.
% Counting starts at 1.
read_header_line(File) :-
    read_line_to_codes(File, Header),
    phrase( header_line(HRN, Samples), Header ),
    abolish( ge_header_name/1 ),
    dynamic( ge_header_name/1 ),
    assertz( ge_header_name(HRN) ),
    abolish( ge_sample_type_participant_pos/4 ),
    dynamic( ge_sample_type_participant_pos/4 ),
    forall(
        nth1(Pos, Samples, Sample),
        (   parse_sample_barcode(Sample, S, T, P),
            assertz( ge_sample_type_participant_pos(S, T, P, Pos) )
        )
    ).

%! read_expression_profile_lines(+File) is det
% Read all expression profile lines and for each line, add a fact
% to the database with the gene identifier and the row of
% normalized gene expression levels for each patient.
read_expression_profile_lines(File) :-
    abolish( gene_expr/2 ),
    dynamic( gene_expr/2 ),
    nb_setval(ge_line_count, 0),
    repeat,
    (   at_end_of_stream(File)
    ->  nb_delete(ge_line_count),
        !
    ;   read_line_to_codes(File, Codes),
        nb_getval(ge_line_count, N),
        (   N mod 1000 =:= 0
        ->  format('.'), flush_output
        ;   true
        ),
        succ(N, N1),
        nb_setval(ge_line_count, N1),
        phrase( expr_line(Gene, Expr), Codes ),
        assertz( gene_expr(Gene, Expr) ),
        fail
    ).

%! parse_sample_barcode(+Sample, -S, -T, -P) is det
% Take a string with a barcode _Sample_, parse it, determining
% its type _T_ and contents, extract the participant information
% _P_ from it, and convert the Sample name to the atom _S_.
parse_sample_barcode(Sample, S, T, P) :-
    phrase( barcode(T, C), Sample ),
    extract_barcode(participant, T, C, P),
    atom_codes(S, Sample).

% DCGs for parsing the gene expression file.
% ------------------------------------------
%
% These are the simple rules for parsing the contents of the two
% different gene expression profile lines. A header line consists
% of a row name _HRN_ and the list of the sample names as strings
% _Ss_, delimited by tabs.
header_line(HRN, Ss) -->
    string_upto_char(0'\t, HeaderRowName),
    {   atom_codes(HRN, HeaderRowName)
    },
    samples(Ss).

% Reading each sample.
samples([S|Ss]) -->
    "\t",
    !,
    string_upto_char(0'\t, S),
    samples(Ss).
samples([]) -->
    [].

% Each expression values line starts with a gene name _Gene_ and
% a list of expression values _Expr_.
expr_line(Gene, Expr) -->
    string_upto_char(0'\t, S),
    exprs(E),
    {   atom_codes(Gene, S),
        Expr =.. [expr|E]
    }.

% When the expression value is missing, there is an "NA" in the
% provided file
exprs([na|Ns]) -->
    "\tNA",
    !,
    exprs(Ns).
% Otherwise, it must be a valid number.
exprs([N|Ns]) -->
    "\t",
    string_upto_char(0'\t, S),
    {   number_codes(N, S)
    }, !,
    exprs(Ns).
exprs([]) -->
    [].

% This reads a string up to a char, much like string_without//2.
% It, however, takes the character _Char_ as an argument and does
% not use memberchk/2 for the comparison.
string_upto_char(Char, [C|T]) -->
    [C],
    {   C \== Char
    }, !,
    string_upto_char(Char, T).
string_upto_char(_Char, []) --> [].


