package satellite

import (
	. "github.com/onsi/ginkgo"
	. "github.com/onsi/gomega"
)

var _ = Describe("spacetrack tests", func() {
	type args struct {
		a spacetrackArgs
	}
	tests := []struct {
		name string
		args args
		want string
	}{
		{
			name: "empty",
			want: "//emptyresult/show",
		},
		{
			name: "some stuff",
			args: args{
				a: spacetrackArgs{
					base:  "base",
					class: "class",
				},
			},
			want: "base/class/emptyresult/show",
		},
		{
			name: "withfilters",
			args: args{a: spacetrackArgs{
				base:  "base",
				class: "class",
				filters: filters{
					{
						filterType: "type",
						operator:   "operator",
						value:      "value",
					},
				},
			}},
			want: "base/class/type/operatorvalue/emptyresult/show",
		},
		{
			name: "multiple filters",
			args: args{a: spacetrackArgs{
				base:  "base",
				class: "class",
				filters: filters{
					{
						filterType: "type",
						operator:   lt,
						value:      "value",
					},
					{
						filterType: "type2",
						operator:   "operator2",
						value:      "value2",
					},
				},
			}},
			want: "base/class/type/<value/type2/operator2value2/emptyresult/show",
		},
		{
			name: "order",
			args: args{a: spacetrackArgs{
				orderByField: "something",
			}},
			want: "//orderby/something asc/emptyresult/show",
		},
		{
			name: "order specified",
			args: args{a: spacetrackArgs{
				orderByField: "something",
				orderByDir:   desc,
			}},
			want: "//orderby/something desc/emptyresult/show",
		},
		{
			name: "limit",
			args: args{a: spacetrackArgs{
				limit: 10,
			}},
			want: "//limit/10/emptyresult/show",
		},
		{
			name: "complete",
			args: args{a: spacetrackArgs{
				base:         baseurl,
				class:        tle,
				orderByField: epoch,
				orderByDir:   desc,
				limit:        1,
				filters: filters{{
					filterType: noradCatID,
					operator:   dq,
					value:      "47966",
				}, {
					filterType: epoch,
					operator:   lt,
					value:      "2021-08-18",
				}},
			}},
			want: "https://www.space-track.org/basicspacedata/query/class/TLE/NORAD_CAT_ID/47966/EPOCH/<2021-08-18/orderby/EPOCH desc/limit/1/emptyresult/show",
		},
	}
	for _, tt := range tests {
		tt := tt
		Describe(tt.name, func() {
			got := buildQuery(tt.args.a)
			It("should match expectations", func() {
				Expect(got).To(Equal(tt.want))
			})
		})
	}
})
