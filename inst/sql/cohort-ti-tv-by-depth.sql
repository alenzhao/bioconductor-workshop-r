# Transition/Transversion Ratio by Depth of Coverage.
SELECT
  (transitions/transversions) AS titv_ratio,
  depth,
FROM (
  SELECT
    SUM(mutation IN ('A->G', 'G->A', 'C->T', 'T->C')) AS transitions,
    SUM(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
                     'A->T', 'T->A', 'C->G', 'G->C')) AS transversions,
    call.DP AS depth,
  FROM (
    SELECT
      CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
      call.DP
    FROM (
      SELECT
        reference_bases,
        GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
        COUNT(alternate_bases) WITHIN RECORD AS num_alts,
        call.DP,
      FROM
        [_THE_TABLE_]
      HAVING
        # Skip 1/2 genotypes and non-SNP variants
        num_alts = 1
        AND reference_bases IN ('A','C','G','T')
        AND alternate_bases IN ('A','C','G','T')
    )
  )
  GROUP BY
    depth
  )
WHERE
  transversions > 0