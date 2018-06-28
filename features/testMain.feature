Feature: Main
  Read input data files, perform computations, write output result files


  Scenario: Provide FASTA files and MutateCompute_Options.txt file
    Given FASTA input files
    When MutateCompute_Options file specifies to perform the Mutate operation
    Then result should be FASTA files for every possible amino acid substitution

    Given FASTA input files
    When MutateCompute_Options file specifies to perform the Agadir operation
    And Agadir/Options.txt specifies to perform Tango
    Then result should be Agadir PSX output files


  Scenario: Provide pdb files and MutateCompute_Options.txt file
    Given pdb input files
    When MutateCompute_Options file specifies to perform the FoldX BuildModel operation
    Then result should be FoldX DDG output files for all specified amino acid substitutions

    Given pdb input files
    When MutateCompute_Options file specifies to perform the FoldX Stability operation
    Then result should be FoldX DDG output files

    Given pdb input files of separate proteins and/or protein chains that form a complex
    When MutateCompute_Options file specifies to perform the FoldX AnalyseComplex operation
    And different proteins of a protein complex are included
    Then result should be FoldX DDG output files for all complexes for which there is a structure
