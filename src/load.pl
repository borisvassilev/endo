/*
This file is part of "endo" Project. It is subject to the license
terms in the LICENSE file found in the top-level directory of
this distribution and at http://opensource.org/licenses/MIT. No
part of "endo" Project, including this file, may be copied,
modified, propagated, or distributed except according to the
terms contained in the LICENSE file.
*/
:- use_module(tcga_data).

load :-
    consult(load),
    consult(tcga_data),
    consult(tcga_barcode).

tests :-
    set_test_options([load(normal),run(manual),silent(false)]),
    run_tests.

data :-
    load_data.

in_data :-
    module(tcga_data).
