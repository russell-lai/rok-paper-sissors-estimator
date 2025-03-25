Estimator for Lattice-based Succinct Arguments
==============================================

The RoK, Paper, SISsors Estimator, or RoK Estimator for short, is a `Sage <http://sagemath.org>`__ module which provides functions for estimating the concrete security and communication costs of lattice-based reduction of knowledge (RoK) protocols. It allows interactively simulate the execution of sequences of RoKs while tracking the parameter changes and costs.

The main purpose of this estimator is allow designers to compare the efficiency and security of lattice-based succinct arguments constructed by different compositions of RoKs known in the literature. 

It has dependencies on the `Lattice Estimator <https://github.com/malb/lattice-estimator>`__ module and some utility functions provided `here <https://github.com/russell-lai/lattice_lib>`__.

Quick Start
-----------

Read the `demo <https://github.com/russell-lai/rok-paper-sissors-estimator/blob/main/rok_estimator_demo.ipynb>`__ to see how to use some basic functions of the module.
 
Status
------

Currently, only RoKs from and to the $\Xi^{\mathsf{lin}}$ relation defined in Section 5 of https://eprint.iacr.org/2024/1972.pdf are covered.

                     
Evolution
---------

This code is evolving, new results are added and bugs are fixed. Hence, estimations from earlier
versions might not match current estimations. This is annoying but unavoidable. We recommend to also
state the commit that was used when referencing this project.

.. warning :: We give no API/interface stability guarantees. We try to be mindful but we may reorganize the code without advance warning.

Bugs
----

Please report bugs through the `GitHub issue tracker <https://github.com/russell-lai/rok-paper-sissors-estimator/issues>`__.

Contributions
-------------

At present, this estimator is maintained by Russell W. F. Lai. Contributors are:

- `Michael Klooß <https://github.com/mklss>`__
- `Russell W. F. Lai <https://github.com/russell-lai>`__
- `Ngoc Khanh Nguyen <https://github.com/khanhcrypto>`__
- `Michał Osadnik <https://github.com/osdnk>`__


Citing
------

If you use this estimator in your work, please cite

    | Michael Klooß, Russell W. F. Lai, Ngoc Khanh Nguyen and Michał Osadnik. *RoK, Paper, SISsors – Toolkit for Lattice-based Succinct Arguments*.
    | Advances in Cryptology – ASIACRYPT 2024. ASIACRYPT 2024. Lecture Notes in Computer Science, vol 15488. Springer, Singapore. 
    | https://doi.org/10.1007/978-981-96-0935-2_7

A more updated version of the above is available as

    | Cryptology ePrint Archive, Report 2024/1972, 2024. https://eprint.iacr.org/2024/1972

License
-------

The estimator is licensed under the `LGPLv3+ <https://www.gnu.org/licenses/lgpl-3.0.en.html>`__ license.

Acknowledgements
----------------

This project was supported through the Research Council of Finland project No. 358951, the Helsinki Institute for Information Technology (HIIT), and the Protocol Labs RFP-013: Cryptonet network grant.
