copy (
with all_regions as (
  SELECT
   region as region,
   CASE WHEN week < 8 * 4 THEN year
        ELSE year + 1
        END as year,
   CASE WHEN week < 8 * 4 THEN week + (52 - 8 * 4) ELSE week - 8 * 4 END
    as week,
   sum("total_patients") over (partition by region, year)
    as Nmax,
   ilitotal as dS
  FROM '{{ input }}'
  WHERE region = 'Region {{ region }}'
  ORDER BY region, year, week
)

SELECT week, dS, Nmax
FROM all_regions
WHERE year = {{ year }}

) to '{{ output }}' with csv header;
