summary.mtsc = function (m, ...)
{	n = m$nr - 1
	t = f = numeric (n)
	for (i in 1:n)
	{	v = m$cp [[i]]
		t [i] = (m$t [v$v1$b] + m$t[v$v2$a]) / 2
		f [i] = 0
	}
	data.frame (t, f)
}

.mtsc.ftest = function (n, nc, ng, rss.f, rss.r)
	(rss.r - rss.f) / rss.f * (n - nc * ng) / nc

