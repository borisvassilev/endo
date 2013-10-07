:- module(markup_source, [make_docs/2]).

:- use_module(library('dcg/basics')).

markup(Options) :-
    forall(
        member(F, Files),
        make_doc(F, asciidoc, Options)
    ).

make_doc(File, asciidoc, Options) :-
    extract_options(Options, RemainingOptions),
    source_markup(user_input, user_output, RemainingOptions),
        
extract_options(Options, Options). % for now!


