# RNG-tests
Testing some RNGs with statistical criteria

Tree statistical criteria with three trust levels for each were used:
1) Uniformity. Wether the bytes of RNG-generated sequence are uniformly distributed.
2) Independency. Wether each byte of sequence is independent of preceding byte.
3) Homogenity. Wether disribution qualities are the same for sufficiently large subsequence of generated sequence of bytes.

Generators tested:

1) Standard C random function.
2) Lehmer (low and high bits).
3) LSR-based generators:
  a) L20;
  b) L89;
  c) Geffe generator.
4) Wolfram generator.
5) Librarian generator (based on original text of 'Moby Dick' by Herman Melville).
