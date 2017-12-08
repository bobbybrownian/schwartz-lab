class CreateMethodtypes < ActiveRecord::Migration[5.1]
  def change
    create_table :methodtypes do |t|

      t.string :method_type, null: false

    end
  end
end
