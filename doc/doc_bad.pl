:- module(source_to_doc, [make_docs/2]).

:- use_module(library('dcg/basics')).

make_docs(Files, Options) :-
    forall(
        member(F, Files),
        make_doc(F, asciidoc, Options)
    ).

make_doc(File, asciidoc, Options) :-
    extract_options(Options, RemainingOptions),
    phrase_from_file( pl_source(SourceTree), File, RemainingOptions ),
    atom_concat(BaseName, '.pl', File),
    atom_concat(BaseName, '.asciidoc', AsciidocFile),
    write_asciidoc(SourceTree, AsciidocFile).
        
extract_options(Options, Options). % for now!

pl_source(CopyRight, Chunks) -->
    empty_lines,
    copyright(CopyRight),
    chunks(Chunks).

empty_lines -->
    blanks_to_nl, !,
    empty_lines.
empty_lines -->
    [].

copyright(Text) -->
    "/*", blanks_to_nl,
    string(Text),
    "*/", blanks_to_nl, !,
    { format('/*~n~s*/~n', [Text]) }.

chunks([H|T]) -->
    chunk(H),
    chunks(T).
chunks([]) --> [].

chunk(doc(D)) -->
    doc(D), !.
chunk(predicate(P)) -->
    predicate(P), !.
chunk(source(S)) -->
    source(S), !.

doc([H|T]) -->
    doc_line(H), !,
    doc(T).
doc([]) --> [].

doc_line(DL) -->
    "% ", string(DL), "\n", !,
    { format('[doc] ~s~n', [DL]) }.

predicate(Decl, Descr, Source) -->
    predicate_declaration(Decl),
    predicate_description(Descr),
    predicate_source(Source).

predicate_declaration(pred_decl(H, D)) -->
    "%! ",
    string(Head),
    blanks, "is",
    blanks, determinism(D),
    blanks, "\n", !,
    {   atom_codes(H, Head)
    }.

determinism(det) -->
    "det".
determinism(semidet) -->
    "semidet".
determinism(failure) -->
    "failure".
determinism(nondet) -->
    "nondet".
determinism(multi) -->
    "multi".
    
predicate_description(D) -->
    doc(D).

predicate_source(S) -->
    source(S).

source([H|T]) -->
    source_line(H), !,
    { format('[src] ~s~n', [H]) },
    source(T).
source([]) --> [].

source_line("") -->
    blanks_to_nl, !.
source_line([0'%,Next|Rest]) -->
    "%", [Next],
    {   \+ code_type(Next, space)
    },
    string(Rest), "\n", !.
source_line([0'%|S]) -->
    [C],
    {   C \== 0'%
    },
    string(S), "\n", !.

