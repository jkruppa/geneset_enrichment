library(shiny)
library(mvtnorm)
library(scatterplot3d)

shinyServer(function(input, output) {

	output$plots <- renderPlot({
		blue2 = rgb(0, 0, 250, maxColorValue=250, alpha=100)
		green2 = rgb(0, 250, 0, maxColorValue=250, alpha=100)

		d = 10000
		path = input$pathgenes
		degs = input$degs
		penr = input$perpath
		degp = degs/d
		enrp = ((penr/100) * path) /path
		alpha = 0.05
		b = 2
		a = seq(0, 100, by=0.1)
		pbet = pbeta(alpha, b, a)
		afinal = a[max(which(pbet<=degp))]
		a2final = a[max(which(pbet<=enrp))]
		p = qbeta(seq(0, 1, length.out=d), b, afinal)
		if (degs==0) p[1] = 0.06
		sum(p<alpha)
		e = qbeta(seq(0, 1, length.out=path), b, a2final)
		e[path] = e[path-1]

		e2 = rep(NA, path)
		for (i in 1:path) e2[i] = min(which(p>=e[i]))
		e = e2
		e = intersect(e, e)
		col.e = rep(3, d)
		col.e[e] = 4
		
		par(mfrow=c(1, 2))
		hist(p, xlim=c(0, 1), col=8, cex.axis=1.5, cex.lab=1.5, xlab="p-value", ylab="frequency", main="", breaks=50)
		points(c(alpha, alpha), c(-100, 1000), type="l", lwd=2, col=2)
		box()
		axis(3, alpha, expression(alpha), cex.axis=1.5, col.ticks=2, lwd.ticks=2)

		plot(-log10(p), cex.lab=1.5, cex.axis=1.5, xlab="rank(p-value)", ylab="-log10(p-value)", col=col.e, pch=21, cex=1.5, lwd=2)
		points((1:d)[which(col.e==4)], -log10(p[which(col.e==4)]), pch=21, col=4, cex=1.5)
		
		points(c(-100, d+100), c(-log10(alpha), -log10(alpha)), type="l", col=2, lwd=2)
		axis(4, -log10(alpha), expression(alpha), cex.axis=1.5, col.ticks=2, lwd.ticks=2)
		legend("topright", c("true", "false"), title="gene in pathway", cex=1.5, col=c(4, 3), pch=15, ncol=2)
	})

	output$plots2 <- renderPlot({
		blue2 = rgb(0, 0, 250, maxColorValue=250, alpha=50)
		green2 = rgb(0, 250, 0, maxColorValue=250, alpha=50)

		d = 1000
		path = input$pathgenes
		degs = input$degs
		penr = input$perpath
		degp = degs/d
		enrp = ((penr/100) * path) /path
		alpha = 0.05
		b = 2
		a = seq(0, 100, by=0.1)
		pbet = pbeta(alpha, b, a)
		afinal = a[max(which(pbet<=degp))]
		a2final = a[max(which(pbet<=enrp))]
		p = qbeta(seq(0, 1, length.out=d), b, afinal)
		if (degs==0) p[1] = 0.06
		sum(p<alpha)
		e = qbeta(seq(0, 1, length.out=path), b, a2final)
		e[path] = e[path-1]

		e2 = rep(NA, path)
		for (i in 1:path) e2[i] = min(which(p>=e[i]))
		e = e2
		e = intersect(e, e)
		col.e = rep(green2, d)
		col.e[e] = blue2
		col.e2 = col.e
		col.e2[col.e==green2] = 0

		par(mfrow=c(1, 2))

		W = wilcox.test(p ~ col.e, alternative="less")
		F = fisher.test(table(p<alpha, col.e))
		
		boxplot(p ~ col.e, cex.lab=1.5, cex.axis=1.5, xlab="gene in pathway", ylab="p-value", names=c("true", "false"), border=c(4, 3), col=c(blue2, green2), main=paste("MWU-test: p =", round(W$p.value, 8)), lwd=2)

		barplot(t(100 * prop.table(table(p<alpha, col.e), 1)), cex.lab=1.5, cex.axis=1.5, xlab="gene significant", ylab="frequency (%)", ylim=c(0, 100), axes=FALSE, col=c(4, 3), cex.names=1.5, names=c("false", "true"), main=paste("Fisher-test: p =", round(F$p.value, 8)))
		axis(2, c(0, 25, 50, 75, 100), c(0, 25, 50, 75, 100), cex.axis=1.5)
	})


})
