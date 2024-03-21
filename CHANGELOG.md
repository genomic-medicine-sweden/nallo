# fellen31/skierfe: Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

<!-- insertion marker -->
<!-- ## [0.1.0](https://github.com/fellen31/skierfe/releases/tag/0.1.0) - 2024-03-21 --> 

### Added

Added test data and test profile [#33](https://github.com/genomic-medicine-sweden/skierfe/pull/33)

---

- Add test profile ([c8922f4](https://github.com/fellen31/skierfe/commit/c8922f415733dbb24944b0461b5160b7538eee87)).
- Add required files and test_data ([ed5b9fa](https://github.com/fellen31/skierfe/commit/ed5b9fa5d6b6ba38df32c4a09d0265ae37d9bdff)).
- add hg38-specific asset files ([2975526](https://github.com/fellen31/skierfe/commit/2975526846c0cac94d97881babd03a30d7a7e433)).
- Add small test data ([7f1a7be](https://github.com/fellen31/skierfe/commit/7f1a7be1eced54fe9fd539674f4a8d8d1f7dedc2)).
- adding codeowners ([36c640b](https://github.com/fellen31/skierfe/commit/36c640b8c1b28d5b80e0f85d000033d3d742f362)).
- Add expected CN-files (for XX/XY), MAF-vcf and option to use an exclude regions file. Also update HiFiCNV version to v.0.1.7 ([90fa4cd](https://github.com/fellen31/skierfe/commit/90fa4cdf424a629c1ef0891826a702400eaa2022)).
- Add grouptuple size so downstream wf's don't have to wait for all samples to complete ([ce696bb](https://github.com/fellen31/skierfe/commit/ce696bb08e22920e79921c21fb400c1cad1c6fc5)).
- Add 'hack' to be able to run ONT-files in TRGT ([00ee540](https://github.com/fellen31/skierfe/commit/00ee5402c97735da51213dd554c2bfe2306a5a26)).
- Add VEP annotation ([261a118](https://github.com/fellen31/skierfe/commit/261a118221962295ca2a858c0da84691fa61f270)).
- Add missing file column to snfs schema ([6cf6b32](https://github.com/fellen31/skierfe/commit/6cf6b3239d896d2bf3f335ef1312c1e79989072a)).
- add config ([23b34b6](https://github.com/fellen31/skierfe/commit/23b34b6eecedf44ce15946cd3e5ea451bffb0473)).
- Add hificnv to workflow ([f4ace8f](https://github.com/fellen31/skierfe/commit/f4ace8f5c5ba3c22a56834efa84ad4c120c1e7be)).
- Add skip for QC ([9910f8c](https://github.com/fellen31/skierfe/commit/9910f8cd76109945e6136b1fb7399450f8bb9036)).
- Add skip for phasing and remove more PED references ([b1277da](https://github.com/fellen31/skierfe/commit/b1277da95b2b6ba040e9f0e7c563a7159e1ba9a9)).
- Add samplesheet, extra_gcvf and extra_snfs validation ([4a8a4c2](https://github.com/fellen31/skierfe/commit/4a8a4c2585b24b5556221b0e46db0387bb763de4)).
- add missing bcftools modules ([0d0cb73](https://github.com/fellen31/skierfe/commit/0d0cb73e490ec8f2810d13d1cbe9220aa48d8141)).
- Add hiphase ([d4ff737](https://github.com/fellen31/skierfe/commit/d4ff737d37dadce7dba687b20eb320b583f75849)).
- Add back --secondary=no to minimap2 align ([178aa08](https://github.com/fellen31/skierfe/commit/178aa089896828d8f9cfb8acafdccf08fa0bc972)).
- Added TRGT for pacbio-data ([3edfdeb](https://github.com/fellen31/skierfe/commit/3edfdeb8fdb89740b04166ee1d43700ad1a4746a)).
- Add fastq and bam stats ([e802266](https://github.com/fellen31/skierfe/commit/e8022662df27ff8c2161120f20dad65a767b831b)).
- Add methylation calls to bedmethyl ([448f9e9](https://github.com/fellen31/skierfe/commit/448f9e9deddef5fbe27b614ae403a198422af30e)).
- add mosdepth ([54e6ea5](https://github.com/fellen31/skierfe/commit/54e6ea5d11aa89116c5a2f7fa038676d9f05e773)).
- Add tandem-regions and reference back to sniffles + skip snv-calling ([c40212c](https://github.com/fellen31/skierfe/commit/c40212c14f8caf02aad909f45cda7a7b1853c28b)).
- Add dipcall as process ([6471dbc](https://github.com/fellen31/skierfe/commit/6471dbc464ef42f83442c8d5d310b83eb9b25732)).
- Add deeptrio, GPU support and fastqc (should get better at keeping features separate) ([64d65ce](https://github.com/fellen31/skierfe/commit/64d65ceb75fb1163aef49664eab35c09931ec62d)).
- Add missing diff..linting fails because of the container URL, not sure how to fix? ([5f55817](https://github.com/fellen31/skierfe/commit/5f55817de626c6278a136d46a172ff383e4dc9da)).
- Added dipcall ([8a8d96e](https://github.com/fellen31/skierfe/commit/8a8d96e3b4411b1d5efaf25c08bd7fc0df6108cb)).
- Added DeepTrio support (incl. GPU for both DV and DT), trios extractedfrom PED-file ([7c5a20a](https://github.com/fellen31/skierfe/commit/7c5a20a26e0011fc4f6015a98f63abb3ea62a31a)).
- Add initial sniffles module ([abc5d14](https://github.com/fellen31/skierfe/commit/abc5d14644636bc48394048ff0932b796c31b7c6)).

### Fixed

- Fix actually skipping CNV workflow is specifying so ([d14a243](https://github.com/fellen31/skierfe/commit/d14a2433513b8fa59e9736f855bda7fd87ea3928)).
- Fix sniffles args bug ([ff5814f](https://github.com/fellen31/skierfe/commit/ff5814ff04fbe1de7f695925f23d3c531ff934a5)).
- Fix mosdepth running without BED ([9c1c5ff](https://github.com/fellen31/skierfe/commit/9c1c5ffb0c73b974316daa6e41a97a9e4102e53e)).
- Fix most process outdirs ([e220943](https://github.com/fellen31/skierfe/commit/e220943675ffe3feef04e7d731d32f53160a9a68)).
- Fix anno output file ([f74269c](https://github.com/fellen31/skierfe/commit/f74269c0f763a9eec4220baa9b2876536e78456e)).
- Fix anno overwriting input with output ([e7d7467](https://github.com/fellen31/skierfe/commit/e7d74678e16a901334a99aa3e622f7ac66a03af2)).
- fix leftover flag in pileup ([5da5dfc](https://github.com/fellen31/skierfe/commit/5da5dfc34b274c4c17f0c5618ee14e0b27a72630)).
- fix methylation bug ([b146e3e](https://github.com/fellen31/skierfe/commit/b146e3e0662c987e6a6ae6eec02b530e46bb14c9)).
- Fix minimap2 index always using default params ([057120d](https://github.com/fellen31/skierfe/commit/057120d67c7dfea214febc471585931555ad8c27)).
- Fix trgt container adress ([38a4d0d](https://github.com/fellen31/skierfe/commit/38a4d0d6cbcd2afbb2f546de729cab975ae4d07c)).
- Fix dipcall sample name in vcf ([85ed680](https://github.com/fellen31/skierfe/commit/85ed680b45757e37046ca5e30c513c5428ca6c1f)).
- Fix cramino phased and unphased ([5a56ef1](https://github.com/fellen31/skierfe/commit/5a56ef10dbbd30b3e2e7c06cbbb47feeb0bd3980)).
- Fix samtools index publishdirs ([a8d45fc](https://github.com/fellen31/skierfe/commit/a8d45fcc02b86096c3d32eaf2b9b514d1ecfb1c2)).
- Fix repeat versions ([571cd30](https://github.com/fellen31/skierfe/commit/571cd30a09f90c5d2dd15c28685221c8e0564ada)).
- fix dipcall not working ([da1e84d](https://github.com/fellen31/skierfe/commit/da1e84d493334c50eb8d90a7f4dd28bbc44df230)).
- Fix methylation sample mix ([8e7a327](https://github.com/fellen31/skierfe/commit/8e7a32756c5569ca2de0abec444b56250b56c5d5)).
- Fix mosdepth not running for all samples ([50a9068](https://github.com/fellen31/skierfe/commit/50a9068e5197b826ae17c5d3bbfaf06504f1cda3)).
- fix dipcall bug ([0789a7d](https://github.com/fellen31/skierfe/commit/0789a7d65e03f2de466c12da8daa1c06fc083d9a)).
- Fix dipcall module? ([46d1c39](https://github.com/fellen31/skierfe/commit/46d1c396ed993f7dca788040fef80e70f11d2910)).
- Fix non-trio hifiasm and add more rememberable skips ([116c0fa](https://github.com/fellen31/skierfe/commit/116c0fa817a9fb2f15fd5f17099315ced67f40fc)).
- fix color ([18a1ba1](https://github.com/fellen31/skierfe/commit/18a1ba187e1ff62f52fbead01cfddac4a926692c)).
- fix minimap2 singularity version ([ba9b67f](https://github.com/fellen31/skierfe/commit/ba9b67fb312a6d2fc09b26c5091833881aca2480)).
- Fix typo ([bd665c3](https://github.com/fellen31/skierfe/commit/bd665c3e9d44e0c94ce0ce920130eaced33bb01d)).

### Changed

- Change to list of dbs ([65a6221](https://github.com/fellen31/skierfe/commit/65a6221bd05e5a6b6e228e15b796b024c8526ffd)).
- Change local GLNexus to nf-core module ([f59e1db](https://github.com/fellen31/skierfe/commit/f59e1db1a01a2bbec386b24d14416043c3f4a383)).
- Change hifiasm to nf-core module, but make it output lowQ beds as well ([84a0e4f](https://github.com/fellen31/skierfe/commit/84a0e4f4394493ff083480741ee792da73c57022)).

### Removed

- Remove nf-core reference ([ada8b25](https://github.com/fellen31/skierfe/commit/ada8b2570240db6cd8e4b35c987f66a839577896)).
- Remove index split alignments & update README ([d841446](https://github.com/fellen31/skierfe/commit/d841446c2cfe36832029b69ae472a81aecd980d5)).
- remove TRVZ because it will create too many processes ([ae8ba8d](https://github.com/fellen31/skierfe/commit/ae8ba8dc821af127599eac326d351632fc6ac9bf)).
- Remove references to non-existing slack ([4e31528](https://github.com/fellen31/skierfe/commit/4e3152847ae11ec3ab12d486a91183af47407906)).
- remove nf-core images and add pipeline summary ([144599e](https://github.com/fellen31/skierfe/commit/144599e198e4075b26ce9ada151a205c55c6ffcf)).
- Remove references to nf-core ([0fbd578](https://github.com/fellen31/skierfe/commit/0fbd578b3aff70258f993e103b1ce79ef442d199)).

### Other

- update path ([0fe838b](https://github.com/fellen31/skierfe/commit/0fe838b45526ecb267ebc4f494aa2aa8805c14b0)).
- Update split_bed_chunks.py ([08cb5ea](https://github.com/fellen31/skierfe/commit/08cb5eaa872e97dacc50485db3ea31fe58b3e868)).
- Update split_bed_chunks.py ([910efca](https://github.com/fellen31/skierfe/commit/910efca54e08b5a91ddba24a210182fa94ee57cc)).
- Update split_bed_chunks.py ([f3593a6](https://github.com/fellen31/skierfe/commit/f3593a60ce3d05496c0c6c639ba57525e8e9722c)).
- Update qc_aligned_reads.nf ([4849d0c](https://github.com/fellen31/skierfe/commit/4849d0cb425304c693c58a574268f1294a6ff69e)).
- Update qc_aligned_reads.nf ([52e9f53](https://github.com/fellen31/skierfe/commit/52e9f53164f2609fd8a306cff6b428c355f0e2c9)).
- Update short_variant_calling.config ([153910f](https://github.com/fellen31/skierfe/commit/153910f96a781e5e4781d2750fa920722fb6ddbf)).
- Update samtools_cat_sort_index.nf ([c048a4e](https://github.com/fellen31/skierfe/commit/c048a4ec7310ab956479a52464aa520b79656609)).
- update hificnv config ([93c9bf8](https://github.com/fellen31/skierfe/commit/93c9bf851a7270d1262b9189d8b6efc24863e13f)).
- Split module configs ([5e1acc6](https://github.com/fellen31/skierfe/commit/5e1acc649eae071682ff40b322cd9ed10668893d)).
- Update docs with pipeline parameters ([6549aa7](https://github.com/fellen31/skierfe/commit/6549aa7fa299109817b142b97f143a174fab51bc)).
- forgot wf ([c510c61](https://github.com/fellen31/skierfe/commit/c510c61b01c535f4a07c98e7b83c5102e336445c)).
- Grab correct one ([c140b55](https://github.com/fellen31/skierfe/commit/c140b5596047173233fbb01733477eb222fdeb67)).
- Singularity image does not contain jemalloc ([4a97504](https://github.com/fellen31/skierfe/commit/4a975047f036a3a46851eca9b31c86f439a486d3)).
- Split alignment for faster processing ([f683091](https://github.com/fellen31/skierfe/commit/f68309102ae6cd4e5d380432c82eae7becb73fda)).
- Change to list of dbs ([65a6221](https://github.com/fellen31/skierfe/commit/65a6221bd05e5a6b6e228e15b796b024c8526ffd)).
- Annotate variants with one database ([6392555](https://github.com/fellen31/skierfe/commit/6392555f6cdaaf21c7f4f1d798e89f447c2caede)).
- Simplify sniffles, and put back dipcall par ([b80d063](https://github.com/fellen31/skierfe/commit/b80d0632a396563b7a904923adffe276b43546b3)).
- update hifiasm ([d78bd1c](https://github.com/fellen31/skierfe/commit/d78bd1c20fde615baf0f673199860b8f640d8c9a)).
- Improve whatshap phasing ([cab48e6](https://github.com/fellen31/skierfe/commit/cab48e629352d104f0fcf64372046df6b1bbd078)).
- minimap index needs to be built separate for dipcall as well ([28fcfdd](https://github.com/fellen31/skierfe/commit/28fcfdd2ed8b635c068c1a47a3f6674c646905cd)).
- Make tandem_repeats optional, since it might not be that great? ([037854b](https://github.com/fellen31/skierfe/commit/037854bffd16e16e544d588df072241e5ce93ce0)).
- Change local GLNexus to nf-core module ([f59e1db](https://github.com/fellen31/skierfe/commit/f59e1db1a01a2bbec386b24d14416043c3f4a383)).
- Change hifiasm to nf-core module, but make it output lowQ beds as well ([84a0e4f](https://github.com/fellen31/skierfe/commit/84a0e4f4394493ff083480741ee792da73c57022)).
- Update dipcall to use PED-file in non-trio mode ([d657b5f](https://github.com/fellen31/skierfe/commit/d657b5f52cc2b37372f3fc1cb6308012b0d2038d)).
- nf-core download would not pull docker image when the container is specified inside single quotes ([990bd79](https://github.com/fellen31/skierfe/commit/990bd79daf4bbb26bb306bac2db0d1bc82277306)).
- See if this way works-3 ([31d0a91](https://github.com/fellen31/skierfe/commit/31d0a913335af96c163f6f1f556b41959c5cb9ac)).
- See if this way works-2 ([16b984d](https://github.com/fellen31/skierfe/commit/16b984d0813acd1a4e284e5073bc3f075534af82)).
- See if this way works ([ead960d](https://github.com/fellen31/skierfe/commit/ead960d053b38829bc44daa0377d65335c597dc2)).
- Update genome_assembly.nf ([c8bfef3](https://github.com/fellen31/skierfe/commit/c8bfef335ff73b5a788b7e47dfbb62bcbec67738)).
- move yak-versions into trio workflow ([10b9121](https://github.com/fellen31/skierfe/commit/10b91218fdb3b4b2780f9c37813483340978473f)).
- Update docs ([31b4ee6](https://github.com/fellen31/skierfe/commit/31b4ee6639b4953332e3b9132c19aa12c89b0d6f)).
- Prettier ([d40c289](https://github.com/fellen31/skierfe/commit/d40c2892477ef7d54796b66c37ae01faa386d3f3)).
- Update modules after TEMPLATE upgrade ([7c2c850](https://github.com/fellen31/skierfe/commit/7c2c850164cf1bf68630528323dc49dcdac93b76)).
- Merge updated TEMPLATE ([dd01545](https://github.com/fellen31/skierfe/commit/dd015455fb30c394d93884568afc4ebb7a5b1df2)).
- Might have messed up git stash ([7c81129](https://github.com/fellen31/skierfe/commit/7c811296d24bc036aa1f08f423f88db230eb6673)).
- Small clean up. ([aca23f3](https://github.com/fellen31/skierfe/commit/aca23f3518198bb65c5a16484a28bfe419d60e1f)).
- linting working for all but samtools/faidx ([238963d](https://github.com/fellen31/skierfe/commit/238963d89c9de2383610f93e072f07782a01d754)).
- Update align.nf ([dc99da2](https://github.com/fellen31/skierfe/commit/dc99da26434f8fc62cbe4f9355adc17c41faccc1)).
- nf-download not working.. ([8d8e9aa](https://github.com/fellen31/skierfe/commit/8d8e9aa1a5cc380485d7a562cd6fa3d567c933ed)).
- nf-download not working.. ([c4c99f4](https://github.com/fellen31/skierfe/commit/c4c99f40b735fbd277a331a71def9ce44779eb7e)).
- Simplify sniffles, and put back dipcall par ([c55f25b](https://github.com/fellen31/skierfe/commit/c55f25b1fce186e1dbd71cd1bd7071dd83ea852a)).
- update hifiasm ([9fde2d4](https://github.com/fellen31/skierfe/commit/9fde2d47eae4e9aef6850007e8f7d8d96b13302c)).
- Improve whatshap phasing ([852c550](https://github.com/fellen31/skierfe/commit/852c550c91855702e13c688b5df17fb5c6b72555)).
- minimap index needs to be built separate for dipcall as well ([0130ee0](https://github.com/fellen31/skierfe/commit/0130ee0dd7e470a1788f2ddb9be0b5b068fd92b8)).
- Make tandem_repeats optional, since it might not be that great? ([28f743d](https://github.com/fellen31/skierfe/commit/28f743ded7f5b93cfd5d44efeda74cde7704d9b7)).
- Change local GLNexus to nf-core module ([2a5f011](https://github.com/fellen31/skierfe/commit/2a5f0117a416af04ac51e52350bcf132147c486f)).
- Change hifiasm to nf-core module, but make it output lowQ beds as well ([2c717d9](https://github.com/fellen31/skierfe/commit/2c717d9c643fbd09cb351e8f6223e94d3286d077)).
- Update dipcall to use PED-file in non-trio mode ([ec3a04d](https://github.com/fellen31/skierfe/commit/ec3a04deada8024862c1165adc6bed4ef83412d0)).
- nf-core download would not pull docker image when the container is specified inside single quotes ([2d023be](https://github.com/fellen31/skierfe/commit/2d023be0526585b34c69742aa20e1837ed2a9733)).
- See if this way works-3 ([5f8a659](https://github.com/fellen31/skierfe/commit/5f8a659c8a0b4aae302adefcefa151033cfc652b)).
- See if this way works-2 ([1a7f47a](https://github.com/fellen31/skierfe/commit/1a7f47a89ee9504556abe466a5c9e98033651154)).
- See if this way works ([595e1cc](https://github.com/fellen31/skierfe/commit/595e1cc8e9660bb118f1cd7ca16a169fb9c48ac7)).
- Update genome_assembly.nf ([27bb008](https://github.com/fellen31/skierfe/commit/27bb008a60b34009f2d6c680a93b0109fad757b2)).
- move yak-versions into trio workflow ([52a2e26](https://github.com/fellen31/skierfe/commit/52a2e261ff710acf666691615ff18ac1594ca862)).
- Update docs ([5ac6516](https://github.com/fellen31/skierfe/commit/5ac6516755f8e8b579c58f13216e22c3d8a338e1)).
- Prettier ([7a0dbd1](https://github.com/fellen31/skierfe/commit/7a0dbd1ac88ed60972a91f2704728f0234811788)).
- Update modules after TEMPLATE upgrade ([5286876](https://github.com/fellen31/skierfe/commit/5286876b53a78c265e95f289574303571af6fefd)).
- Merge updated TEMPLATE ([06a68c0](https://github.com/fellen31/skierfe/commit/06a68c033c745564acb86892aaa26b8b2ec4c288)).
- Might have messed up git stash ([78fe1f5](https://github.com/fellen31/skierfe/commit/78fe1f59d131ba745fcc7fdd4cfc13fd9fc8608a)).
- Small clean up. ([2fb2fea](https://github.com/fellen31/skierfe/commit/2fb2fea11b71ef231c230b17f0dc6d8f2434f201)).
- linting working for all but samtools/faidx ([2fd4fd2](https://github.com/fellen31/skierfe/commit/2fd4fd20e328543b21f5ee292babfc2d57634bb4)).
- Update align.nf ([db6c8e3](https://github.com/fellen31/skierfe/commit/db6c8e3c706ae4ddf5a3a08d200fbd7579de0fa8)).
- More accurate GLNexus process label ([70fc63a](https://github.com/fellen31/skierfe/commit/70fc63a55196c229053297d6be106ebdeb7c9361)).
- More accurate Deepvariant process label ([8be2d12](https://github.com/fellen31/skierfe/commit/8be2d12a8d83980c7a1cc85795c46d899d7fcf3d)).
- Small clean up. ([f0e5ebb](https://github.com/fellen31/skierfe/commit/f0e5ebb51f6df2e5de3b5ac00b69728b56a339de)).
- Modified samplesheet header names for consistent naming ([752f115](https://github.com/fellen31/skierfe/commit/752f11548da6fc96fa60d2f6403760859827c58b)).
- Think working with pbmm2 alignment, uploading to check on uppmax ([f162402](https://github.com/fellen31/skierfe/commit/f16240262396d439467a48aeaf0619d7d898cbaf)).
- initial template build from nf-core/tools, version 2.7.2 ([e84b30a](https://github.com/fellen31/skierfe/commit/e84b30a8b73d246e14afdebf4076ea5d38346197)).
