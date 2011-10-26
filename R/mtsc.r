mtsc = function (t, mts, fit=TRUE, standardise=TRUE)
{	if (!is.matrix (mts) ) mts = as.matrix (mts)
	if (standardise) mts = .mtsc.standardise (mts)
	m = extend (new.env (), "mtsc", t, mts, clusters=NA, cp=NA,
		nr=nrow (mts), nc=ncol (mts) ) 
	if (fit) .fit.mtsc (m)
	m
}

.fit.mtsc = function (m, ...)
{	chain1 = .mtsc.chain.primary (m)
	chain2 = .mtsc.chain.secondary (m, chain1)
	distance = objref (numeric (m$nr - 1) )
	for (i in 1:(m$nr - 1) )
		distance [i] = chain2 [[i]]$rss

	m$cp = vector ("list", m$nr - 1)
	while (length (chain1) > 1)
		.mtsc.join (m, chain1, chain2, distance)

	m$clusters = .mtsctree (chain1 [[1]])

	invisible (m)
}

.mtsc.join = function (m, chain1, chain2, distance)
{	i = which.min (distance () )

	v1 = chain1 [[i]]
	v2 = chain1 [[i + 1]]
	u = chain2 [[i]]
	chain1 [[i]] = .complexnode (v1, v2, v1$a, v2$b, u$mean, u$rss)
	chain1 [[i + 1]] = NULL
	chain2 [[i]] = NULL
	distance$obj = distance [-i]

	m$cp [[length (chain1)]] = chain1 [[i]]

	if (i > 1)
	{	chain2 [[i - 1]] = .mtsc.prec (m, chain1 [[i - 1]], chain1 [[i]])
		distance [i - 1] = chain2 [[i - 1]]$rss -
			chain1 [[i - 1]]$rss - chain1 [[i]]$rss

	}
	if (i < length (chain1) )
	{	chain2 [[i]] = .mtsc.prec (m, chain1 [[i]], chain1 [[i + 1]])
		distance [i] = chain2 [[i]]$rss -
			chain1 [[i]]$rss - chain1 [[i + 1]]$rss
	}
}

.mtsc.chain.primary = function (m)
{	chain = vector ("list", m$nr)
	for (i in 1:m$nr) chain [[i]] = .simplenode (i, mean=m$mts [i,])
	objref (chain)	
}

.mtsc.chain.secondary = function (m, chain1)
{	chain = vector ("list", m$nr - 1)
	for (i in 1:(m$nr - 1) )
	{	mean = (chain1 [[i]]$mean + chain1 [[i + 1]]$mean) / 2
		rss = .mtsc.rss (m$mts [c (i, i + 1),], mean)
		chain [[i]] = list (mean=mean, rss=rss)		
	}
	objref (chain)
}

.mtsctree = function (v)
	extend (v, "mtsctree")

.simplenode = function (a, b=a, mean=0, rss=0)
	implant (new.env (), a, b, mean, rss)

.complexnode = function (v1, v2, a, b, mean=0, rss=0)
	implant (new.env (), v1, v2, a, b, mean, rss)

.mtsc.prec = function (m, v1, v2)
{	n1 = v1$b - v1$a + 1
	n2 = v2$b - v2$a + 1
	a = v1$a
	b = v2$b
	mean = (n1 * v1$mean + n2 * v2$mean) / (n1 + n2)
	rss = .mtsc.rss (m$mts [v1$a:v2$b,], mean)
	list (mean=mean, rss=rss)
}

.mtsc.standardise = function (mts)
{	nc = ncol (mts)
	for (j in 1:nc)
	{	mts [,j] = mts [,j] - mean (mts [,j])
		mts [,j] = mts [,j] / mean (abs (mts [,j]) )
	}
	mts
}

.mtsc.rss = function (x, xbar)
{	y = 0
	if (nrow (x) <= ncol (x) ) for (i in 1:nrow (x) ) y = y + sum ( (x [i,] - xbar)^2)
	else for (i in 1:length (xbar) ) y = y + sum ( (x [,i] - xbar [i])^2)
	y
}





