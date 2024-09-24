# Issues
CDR-Finder may fail to correctly call the correct regions. Choosing the correct set of parameters may fix similar scenarios.

### Edge regions in pericentromeric regions
Regions on the edges may get falsely called.
![](issues/HG00268_haplotype2-0000071.png)
![](issues/HG00358_haplotype2-0000083.png)
![](issues/HG01457_haplotype2-0000112.png)


### Low, uniform average methylation percent
Using the default parameters may not accurately call CDRs. For example, low and uniform average methylation percent may result in false calls.

With a `prom_perc_valley_threshold` of `0.2`.
![](issues/NA19331_haplotype1-0000011_prom0.2.png)


Increasing `prom_perc_valley_threshold` to `0.3` correctly ignores non-prominent CDRs.
![](issues/NA19331_haplotype1-0000011_prom0.3.png)
