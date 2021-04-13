# Dairying, disease and the evolution of lactase persistence in Europe

## Process phenotypes

```sh
# start job
salloc --nodes=1 --cpus-per-task=21 --mem=80G --time=06:00:00 --partition=mrcieu
# use node ramdisk for faster processing
d=$(mktemp -d)
echo "copying pheno file to $d"
cp /mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/16729/2020-11-13/data/data.43017.tab "$d"/
echo "done"
```

## Forest plot

```sh
Rscript assoc.R
```

## Hardy-Weinberg test

```sh
Rscript hwe.R
```

## Spousal analysis

```sh
Rscript spousal.R
```