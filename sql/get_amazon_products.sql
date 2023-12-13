-- Choose a sample of 100 amazon products matching specific criteria
copy (
with reviews as (
  SELECT
    id,
    ts - min(ts) over (partition by id) as ts
  FROM read_csv(
    'Electronics.csv',
    header=false,
    columns={
    'id': 'VARCHAR',
    'reviewer_id': 'VARCHAR',
    'rating': 'BIGINT',
    'ts': 'BIGINT'
  })
)

, valid_products as (
  SELECT
    id,
    count(CASE WHEN ts < 60 * 60 * 24 * 365 * 4 THEN 1 END)
       as reviews_4y,
    max(ts) as last_review_ts
  FROM reviews
  GROUP BY id
  HAVING reviews_4y > 100
     AND last_review_ts > 60 * 60 * 24 * 365 * 4
)

SELECT *
  FROM valid_products
 using sample reservoir(100 rows) repeatable (42)
) to 'amazon-products.csv' with csv header;
