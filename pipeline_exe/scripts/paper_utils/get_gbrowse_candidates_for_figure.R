r = read.delim("aligned.peaks.annotated.capped_root.filtered", sep = ",", header=T)
l = read.delim("aligned.peaks.annotated.capped_leaf.filtered", sep = ",", header=T)

a = merge(r,l, by = "TranscriptID", suffixes=c(".root", ".leaf"))
a$dist = abs(a$ModeLocation.root - a$ModeLocation.leaf)


m = aggregate(a$GeneName.root, list(a$TranscriptID), length)
m = m[m$x < 3,]
a = a[a$TranscriptID %in% m$Group.1,]

p = 400

# both TSS strong and in diff locations
a1 =  a[a$ModeReadCount.root > p & a$ModeReadCount.leaf > p & a$dist > 40 ,]
a1 = a1[,c(1,6,7,8,18:20,26)]

a2 = a[a$ModeReadCount.root > p & a$ModeReadCount.leaf > p & a$dist < 10 & a$Shape.root == "NarrowPeak" ,]

