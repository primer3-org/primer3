DESTDIR ?= build

targets = ntdpal ntthal oligotm primer3_core primer3_masker

base:
	cd src && make

test: base
	cd src && make test

install: base
	mkdir -p $(DESTDIR)/bin
	cp $(addprefix src/, $(targets)) $(DESTDIR)/bin
	cp bin/primer3 $(DESTDIR)/bin
	cp -r share $(DESTDIR)

uninstall:
	rm -f $(addprefix $(DESTDIR)/bin/, $(targets))
	rm -f $(DESTDIR)/bin/primer3
	rm -rf $(DESTDIR)/share/primer3

clean:
	cd src && make clean

