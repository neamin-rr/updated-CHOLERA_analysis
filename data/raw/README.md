# Input data

These files are **not** included in the repository. The surveillance data are held by the
Ethiopian Public Health Institute and are available on reasonable request.

Place here:

| File | Description |
|---|---|
| `cholera_zone_month_2015_2024.csv` | zone, year, month, cases, monthly_population - one row per zone-month (92 x 120 = 11,040 rows) |
| `Adm_2/eth_admbnda_adm2_csa_bofedb_2021.*` | OCHA Ethiopia admin-2 boundaries (shp, shx, dbf, prj) |

The boundary files are public: https://data.humdata.org/dataset/cod-ab-eth

Column names expected by `R/01_stage1_hurdle_inla.R`:

```
zone                 character, admin-2 name
year                 integer, 2015-2024
month                integer, 1-12
cases                integer, reported cholera cases
monthly_population   numeric, zonal population for that month
```
