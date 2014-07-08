var coge = (function (namespace) {
    namespace.Grid = function(container, options, columns) {
        var self = this;

        this.container = $(container),
        this.columns = columns || {},
        this.options = options || {},
        this.options.filterArgs = options.filterArgs || {
            show: 0,
            searchType: 1,
            searchString: '',
        };
        this.options.inlineFilters = options.inlineFilters || true;

        this.options.selectionModel = options.selectionModel ||
        new Slick.RowSelectionModel({
            selectActiveRow: false
        });

        this.dataView = new Slick.Data.DataView(this.options.dataView);
        this.grid = new Slick.Grid(container, this.dataView, this.columns,
            this.options);

        // Grid callbacks
        this.grid.setSelectionModel(this.options.selectionModel);

        this.grid.onCellChange.subscribe(function(e, args) {
        });

        this.grid.onSort.subscribe(function(e, args) {
            //if (e) return console.warn(e);
            //
            var comparator = function(a, b) {
                var index = args.sortCol.field;
                return self.options.comparator(a[index], b[index]);
            }

            self.dataView.sort(comparator, args.sortAsc);
        });

        /// Dataview callbacks
        this.dataView.onRowsChanged.subscribe(function(e, args) {
            self.grid.invalidateRows(args.rows);
            self.grid.render();
        });

        this.dataView.onRowCountChanged.subscribe(function(e, args) {
            self.grid.updateRowCount();
            self.grid.render();
            if (self.container.is(":not(visible)")) {
                self.container.slideDown();
            }
        });

        // Column picker
        this.columnPicker = new Slick.Controls.ColumnPicker(this.columns,
                this.grid, this.options);

        if (this.options.filter) {
            this.dataView.setFilter(this.options.filter);
            this.dataView.setFilterArgs(this.options.filterArgs);
        }

        this.dataView.refresh();
    };

    namespace.Grid.prototype.load = function(data) {
        this.dataView.beginUpdate();
        this.dataView.setItems(data);
        this.dataView.setFilterArgs(this.filterArgs);

        if (this.options.filter) {
            this.dataView.setFilter(this.options.filter);
            this.dataView.setFilterArgs(this.options.filterArgs);
        }

        this.dataView.endUpdate();
        this.dataView.refresh();
        this.filter();
    };

    namespace.Grid.prototype.filter = function() {
        if (!this.options.filter) return;

        this.dataView.refresh();
    }

    return namespace;

})(coge || {});
