Feature: Main
  Read input data files, perform computations, write output result files

  Scenario: Provide FASTA files and MutateCompute_Options.txt file
    Given FASTA files
    When MutateCompute_Options file specifies to perform the Mutate operation
    Then result should be FASTA files for every possible amino acid substitution

    Given FASTA files
    When MutateCompute_Options file specifies to perform the Agadir operation
    And Agadir/Options.txt specifies to perform Tango
    Then result should be Agadir PSX output files

  Scenario: Provide pdb files and MutateCompute_Options.txt file
    Given FASTA files
    When MutateCompute_Options file specifies to perform the Mutate operation
    Then result should be FASTA files for every possible amino acid substitution
