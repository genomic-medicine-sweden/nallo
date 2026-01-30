BEGIN { OFS="\t" }

/^##/ { print; next }

/^#CHROM/ {
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,"sample"
    next
}

{
    split($9, fmt, ":")
    nfmt = length(fmt)

    gt = ""
    dp_sum = 0
    delete ad_sum

    for (i=10; i<=NF; i++) {
        split($i, f, ":")

        for (j=1; j<=nfmt; j++) {
            if (fmt[j] == "GT") {
                gt = gt (gt ? "/" : "") f[j]
            }
            else if (fmt[j] == "DP") {
                dp_sum += f[j]
            }
            else if (fmt[j] == "AD") {
                split(f[j], ad, ",")
                for (k=1; k<=length(ad); k++) {
                    ad_sum[k] += ad[k]
                }
            }
        }
    }

    # rebuild FORMAT field
    out = ""
    for (j=1; j<=nfmt; j++) {
        if (fmt[j] == "GT") {
            out = out (out ? ":" : "") gt
        }
        else if (fmt[j] == "DP") {
            out = out (out ? ":" : "") dp_sum
        }
        else if (fmt[j] == "AD") {
            ads = ""
            for (k=1; k<=length(ad_sum); k++) {
                ads = ads (ads ? "," : "") ad_sum[k]
            }
            out = out (out ? ":" : "") ads
        }
        else {
            # drop unsupported FORMAT fields silently
            out = out (out ? ":" : "") "."
        }
    }

    $10 = out
    NF = 10
    print
}
