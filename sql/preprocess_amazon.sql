copy (
with base as (
  SELECT
    (ts - min(ts) over ( partition by id ))
      // (60 * 60 * 24 * 7) as week
  FROM read_csv(
    'Electronics.csv',
    header=false,
    columns={
    'id': 'VARCHAR',
    'reviewer_id': 'VARCHAR',
    'rating': 'BIGINT',
    'ts': 'BIGINT'
  })
  WHERE id = '{{ id }}'
)

SELECT
  week, count(1) as dS, 1e5 as Nmax
from base
group by week
order by week

) to '{{ output }}' with csv header;
