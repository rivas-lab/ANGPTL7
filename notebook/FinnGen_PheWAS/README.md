# FinnGen PheWAS

The data is downloaded from:
http://results.finngen.fi/variant/1-11193760-C-T

## How to uncrop the svg file?

```
cat locuszoom.svg  | sed -e 's/font-size: 12px/font-size: 14px/g' | sed -e 's/font-size: 11px/font-size: 10px/g' | sed -e 's/554/900/g' | sed -e 's/1478/1600/g' > locuszoom.full.svg
```
