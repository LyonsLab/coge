from cdagline import DagLine
from evalue_adjust import adjust_evalue

def main():
    import sys

    commands = ('plot', 'from_blast', 'adjust')
    if len(sys.argv) < 2 or not sys.argv[1] in commands:
        print "dagtools\n"
        print "available commands are: " 
        for c in commands:
            print "\t", c
        print "\ne.g as '%s plot -h' to see the help for a particular command" \
                % (sys.argv[0])
        sys.exit()

    command = sys.argv[1]
    if command == 'plot':
        import plot_dag
        plot_dag.main(sys.argv[2:])
    elif command == 'from_blast':
        import from_blast
        from_blast.main(sys.argv[2:])
    elif command == "adjust":
        import evalue_adjust
        evalue_adjust.main(sys.argv[2:])


if __name__ == "__main__":
    main()
